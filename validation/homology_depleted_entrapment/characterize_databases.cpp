// Pre-search characterization of the five preregistered target databases.
// Validation-only; reads FASTA and removal manifests, never search outcomes.
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr int NCOND = 5;
constexpr int HLL_P = 18;
constexpr std::size_t HLL_M = 1ULL << HLL_P;
constexpr double STATIC_C = 57.021464;
constexpr double WATER_PLUS_PROTON = 19.018389715;
constexpr std::array<double, 10> MASS_EDGES = {600, 800, 1000, 1200, 1500, 2000,
                                               2500, 3000, 4000, 5000};
const std::array<std::string, NCOND> NAMES = {
    "original", "homology_depleted", "size_control_130363",
    "size_control_155921", "size_control_196613"};

struct Record {
  std::string header;
  std::string sequence;
  bool entrapment = false;
  int entrapment_index = -1;
};

uint64_t mix64(uint64_t value) {
  value ^= value >> 30;
  value *= 0xbf58476d1ce4e5b9ULL;
  value ^= value >> 27;
  value *= 0x94d049bb133111ebULL;
  return value ^ (value >> 31);
}

struct Hll {
  std::vector<uint8_t> registers = std::vector<uint8_t>(HLL_M, 0);

  void add(uint64_t input) {
    uint64_t hash = mix64(input);
    std::size_t index = hash & (HLL_M - 1);
    uint64_t word = hash >> HLL_P;
    int rank = word ? (__builtin_clzll(word) - HLL_P + 1) : (64 - HLL_P + 1);
    registers[index] = std::max(registers[index], static_cast<uint8_t>(rank));
  }

  void merge(const Hll& other) {
    for (std::size_t i = 0; i < HLL_M; ++i) {
      registers[i] = std::max(registers[i], other.registers[i]);
    }
  }

  double estimate() const {
    double sum = 0.0;
    std::size_t zeros = 0;
    for (uint8_t value : registers) {
      sum += std::ldexp(1.0, -static_cast<int>(value));
      zeros += value == 0;
    }
    double m = static_cast<double>(HLL_M);
    double alpha = 0.7213 / (1.0 + 1.079 / m);
    double raw = alpha * m * m / sum;
    if (raw <= 2.5 * m && zeros) {
      return m * std::log(m / static_cast<double>(zeros));
    }
    return raw;
  }
};

struct Stats {
  uint64_t proteins = 0;
  uint64_t residues = 0;
  uint64_t standard_residues = 0;
  std::array<uint64_t, 26> amino{};
  uint64_t fully_tryptic_instances = 0;
  uint64_t searchable_instances = 0;
  uint64_t one_enzymatic_terminus = 0;
  uint64_t two_enzymatic_termini = 0;
  std::array<uint64_t, 64> searchable_length{};
  std::array<uint64_t, MASS_EDGES.size() - 1> searchable_mass{};
  Hll fully_tryptic_unique;
  Hll searchable_unique;

  void merge(const Stats& other) {
    proteins += other.proteins;
    residues += other.residues;
    standard_residues += other.standard_residues;
    fully_tryptic_instances += other.fully_tryptic_instances;
    searchable_instances += other.searchable_instances;
    one_enzymatic_terminus += other.one_enzymatic_terminus;
    two_enzymatic_termini += other.two_enzymatic_termini;
    for (std::size_t i = 0; i < amino.size(); ++i) amino[i] += other.amino[i];
    for (std::size_t i = 0; i < searchable_length.size(); ++i)
      searchable_length[i] += other.searchable_length[i];
    for (std::size_t i = 0; i < searchable_mass.size(); ++i)
      searchable_mass[i] += other.searchable_mass[i];
    fully_tryptic_unique.merge(other.fully_tryptic_unique);
    searchable_unique.merge(other.searchable_unique);
  }
};

bool standard(char residue) {
  static const std::string alphabet = "ACDEFGHIKLMNPQRSTVWY";
  return alphabet.find(residue) != std::string::npos;
}

double residue_mass(char residue) {
  switch (residue) {
    case 'A': return 71.037113805;
    case 'R': return 156.101111050;
    case 'N': return 114.042927470;
    case 'D': return 115.026943065;
    case 'C': return 103.009184505 + STATIC_C;
    case 'E': return 129.042593135;
    case 'Q': return 128.058577540;
    case 'G': return 57.021463735;
    case 'H': return 137.058911875;
    case 'I': case 'L': return 113.084063975;
    case 'K': return 128.094963015;
    case 'M': return 131.040484645;
    case 'F': return 147.068413945;
    case 'P': return 97.052763875;
    case 'S': return 87.032028435;
    case 'T': return 101.047678505;
    case 'W': return 186.079312980;
    case 'Y': return 163.063328575;
    case 'V': return 99.068413945;
    default: return -1.0;
  }
}

int mass_bin(double mass) {
  for (std::size_t i = 0; i + 1 < MASS_EDGES.size(); ++i) {
    if (mass >= MASS_EDGES[i] && mass < MASS_EDGES[i + 1]) return static_cast<int>(i);
  }
  if (mass == MASS_EDGES.back()) return static_cast<int>(MASS_EDGES.size() - 2);
  return -1;
}

std::vector<Record> read_fasta(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open FASTA " + path);
  std::vector<Record> records;
  std::string line, header, sequence;
  int entrapment_index = 0;
  auto emit = [&]() {
    if (header.empty()) return;
    bool entrapment = header.rfind("ENT_", 0) == 0;
    records.push_back({header, sequence, entrapment,
                       entrapment ? entrapment_index++ : -1});
  };
  while (std::getline(input, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty() && line[0] == '>') {
      emit();
      header = line.substr(1);
      sequence.clear();
    } else {
      sequence += line;
    }
  }
  emit();
  return records;
}

std::unordered_set<int> read_removals(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("cannot open removal file " + path);
  std::unordered_set<int> result;
  std::string line;
  std::getline(input, line);
  while (std::getline(input, line)) {
    std::size_t tab = line.find('\t');
    result.insert(std::stoi(line.substr(0, tab)));
  }
  return result;
}

void update_candidate(std::array<Stats, NCOND>& total,
                      std::array<Stats, NCOND>& entrapment,
                      const std::array<bool, NCOND>& active, bool is_entrapment,
                      int start, int end, double mass, bool fully_tryptic,
                      const std::vector<uint64_t>& prefix_hash,
                      const std::vector<uint64_t>& powers) {
  int length = end - start;
  if (length < 1 || length > 63 || mass < 600.0 || mass > 5000.0) return;
  uint64_t hash = prefix_hash[end] - prefix_hash[start] * powers[length];
  hash ^= static_cast<uint64_t>(length) * 0x9e3779b97f4a7c15ULL;
  int mbin = mass_bin(mass);
  for (int condition = 0; condition < NCOND; ++condition) {
    if (!active[condition]) continue;
    Stats* destinations[2] = {&total[condition], is_entrapment ? &entrapment[condition] : nullptr};
    for (Stats* stats : destinations) {
      if (!stats) continue;
      ++stats->searchable_instances;
      ++stats->searchable_length[length];
      if (mbin >= 0) ++stats->searchable_mass[mbin];
      stats->searchable_unique.add(hash);
      if (fully_tryptic) {
        ++stats->fully_tryptic_instances;
        ++stats->two_enzymatic_termini;
        stats->fully_tryptic_unique.add(hash);
      } else {
        ++stats->one_enzymatic_terminus;
      }
    }
  }
}

void characterize_record(const Record& record,
                         const std::array<std::unordered_set<int>, NCOND>& removed,
                         std::array<Stats, NCOND>& total,
                         std::array<Stats, NCOND>& entrapment) {
  std::array<bool, NCOND> active{};
  for (int condition = 0; condition < NCOND; ++condition) {
    active[condition] = !record.entrapment ||
                        removed[condition].find(record.entrapment_index) == removed[condition].end();
    if (!active[condition]) continue;
    ++total[condition].proteins;
    total[condition].residues += record.sequence.size();
    if (record.entrapment) {
      ++entrapment[condition].proteins;
      entrapment[condition].residues += record.sequence.size();
    }
    for (char residue : record.sequence) {
      if (residue >= 'A' && residue <= 'Z') {
        ++total[condition].amino[residue - 'A'];
        if (record.entrapment) ++entrapment[condition].amino[residue - 'A'];
      }
      if (standard(residue)) {
        ++total[condition].standard_residues;
        if (record.entrapment) ++entrapment[condition].standard_residues;
      }
    }
  }

  const std::string& sequence = record.sequence;
  int n = static_cast<int>(sequence.size());
  std::vector<int> boundaries{0};
  for (int i = 0; i < n; ++i) {
    if ((sequence[i] == 'K' || sequence[i] == 'R') &&
        (i + 1 == n || sequence[i + 1] != 'P')) {
      boundaries.push_back(i + 1);
    }
  }
  if (boundaries.back() != n) boundaries.push_back(n);
  std::vector<int> boundary_index(n + 1, -1);
  for (int i = 0; i < static_cast<int>(boundaries.size()); ++i)
    boundary_index[boundaries[i]] = i;

  constexpr uint64_t BASE = 1315423911ULL;
  std::vector<uint64_t> prefix_hash(n + 1, 0), powers(64, 1);
  for (int i = 1; i < 64; ++i) powers[i] = powers[i - 1] * BASE;
  for (int i = 0; i < n; ++i) {
    char residue = sequence[i] == 'I' ? 'L' : sequence[i];
    prefix_hash[i + 1] = prefix_hash[i] * BASE + static_cast<unsigned char>(residue) + 1;
  }

  // Peptides with an enzymatic N terminus. This also emits all fully tryptic
  // peptides, including up to two internal missed cleavages.
  for (int bi = 0; bi + 1 < static_cast<int>(boundaries.size()); ++bi) {
    int start = boundaries[bi];
    int limit = std::min({n, start + 63,
                          boundaries[std::min(bi + 3, static_cast<int>(boundaries.size()) - 1)]});
    double mass = WATER_PLUS_PROTON;
    for (int end = start + 1; end <= limit; ++end) {
      double value = residue_mass(sequence[end - 1]);
      if (value < 0) break;
      mass += value;
      if (mass > 5000.0) break;
      bool full = boundary_index[end] >= 0;
      update_candidate(total, entrapment, active, record.entrapment, start, end,
                       mass, full, prefix_hash, powers);
    }
  }

  // Peptides with only an enzymatic C terminus. Starts that are enzymatic are
  // skipped because the preceding loop already emitted them.
  for (int bj = 1; bj < static_cast<int>(boundaries.size()); ++bj) {
    int end = boundaries[bj];
    int limit = std::max({0, end - 63, boundaries[std::max(0, bj - 3)]});
    double mass = WATER_PLUS_PROTON;
    for (int start = end - 1; start >= limit; --start) {
      double value = residue_mass(sequence[start]);
      if (value < 0) break;
      mass += value;
      if (mass > 5000.0) break;
      if (boundary_index[start] >= 0) continue;
      update_candidate(total, entrapment, active, record.entrapment, start, end,
                       mass, false, prefix_hash, powers);
    }
  }
}

void print_array(std::ostream& out, const auto& values) {
  out << '[';
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) out << ',';
    out << values[i];
  }
  out << ']';
}

void print_stats(std::ostream& out, const Stats& stats) {
  out << "{\"proteins\":" << stats.proteins
      << ",\"residues\":" << stats.residues
      << ",\"standard_residues\":" << stats.standard_residues
      << ",\"fully_tryptic_peptide_instances\":" << stats.fully_tryptic_instances
      << ",\"fully_tryptic_unique_hll\":" << std::llround(stats.fully_tryptic_unique.estimate())
      << ",\"searchable_semi_tryptic_peptide_instances\":" << stats.searchable_instances
      << ",\"searchable_unique_hll\":" << std::llround(stats.searchable_unique.estimate())
      << ",\"one_enzymatic_terminus\":" << stats.one_enzymatic_terminus
      << ",\"two_enzymatic_termini\":" << stats.two_enzymatic_termini
      << ",\"amino_acid_counts_A_to_Z\":";
  print_array(out, stats.amino);
  out << ",\"searchable_length_counts_0_to_63\":";
  print_array(out, stats.searchable_length);
  out << ",\"searchable_mass_edges\":";
  print_array(out, MASS_EDGES);
  out << ",\"searchable_mass_counts\":";
  print_array(out, stats.searchable_mass);
  out << '}';
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 7) {
    std::cerr << "usage: characterize_databases ORIGINAL.fasta HOM.removed.tsv "
                 "C1.removed.tsv C2.removed.tsv C3.removed.tsv OUTPUT.json\n";
    return 2;
  }
  try {
    auto records = read_fasta(argv[1]);
    std::array<std::unordered_set<int>, NCOND> removed;
    for (int condition = 1; condition < NCOND; ++condition)
      removed[condition] = read_removals(argv[condition + 1]);

    int threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif
    std::vector<std::array<Stats, NCOND>> local_total(threads);
    std::vector<std::array<Stats, NCOND>> local_entrapment(threads);
    std::atomic<uint64_t> completed{0};
#pragma omp parallel for schedule(dynamic, 64)
    for (std::size_t index = 0; index < records.size(); ++index) {
      int thread = 0;
#ifdef _OPENMP
      thread = omp_get_thread_num();
#endif
      characterize_record(records[index], removed, local_total[thread], local_entrapment[thread]);
      uint64_t done = ++completed;
      if (done % 50000 == 0) {
#pragma omp critical
        std::cerr << "characterized " << done << '/' << records.size() << " proteins\n";
      }
    }

    std::array<Stats, NCOND> total, entrapment;
    for (int thread = 0; thread < threads; ++thread) {
      for (int condition = 0; condition < NCOND; ++condition) {
        total[condition].merge(local_total[thread][condition]);
        entrapment[condition].merge(local_entrapment[thread][condition]);
      }
    }
    std::ofstream output(argv[6]);
    output << std::setprecision(12)
           << "{\n \"schema_version\":1,\n"
           << " \"searchable_definition\":\"semi-tryptic, 0-2 missed cleavages, length 1-63, MH+ 600-5000, static C+57.021464, no ambiguous residues\",\n"
           << " \"unique_method\":\"HyperLogLog p=18 over I/L-canonical 64-bit rolling hashes; nominal relative SE 0.203%\",\n"
           << " \"conditions\":{\n";
    for (int condition = 0; condition < NCOND; ++condition) {
      if (condition) output << ",\n";
      output << "  \"" << NAMES[condition] << "\":{\"complete_target_database\":";
      print_stats(output, total[condition]);
      output << ",\"entrapment_component\":";
      print_stats(output, entrapment[condition]);
      output << '}';
    }
    output << "\n }}\n";
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
