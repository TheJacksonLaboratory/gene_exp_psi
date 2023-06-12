// For more information, please refer to section "Enriched motif testing" in the
// Methods.
// Build the executable using the command: g++ -pthread -o permute permute.cc
// Retrieve and extract motif locations from:
// genecascade.org/downloads/gene_exp_psi/motif_locations.zip
// Available values for repeat_id are as follows: 0 (CPEs), 1 (TFFMs), 2 (RBPs)

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/wait.h>
#include <unistd.h>
#endif

using namespace std;

void it_step(int i_begin, int i_end, int thread_no, unsigned my_seed);
string group_str(int i);
bool fexists(string fname);
string run_time_str(int no);
void pnt_top_summary(int cnt_tests, double adj_p, string run_time_strs[]);
void pnt_res(int i, int pos, int thread_no);

unsigned no_seqs[3] = {0, 0, 0};

ofstream outfile;
#ifndef _DEBUG
int iterations = 10000;
int no_threads = 7;
#else
int iterations = 10;
int no_threads = 1;
#endif

bool repeat_run = true;

// p-values equal to or less than the threshold are marked with an asterisk (*)
// for emphasis.
const double p_threshold = 0.05;
// Tests where the relative or absolute difference between base 0 and base UHP
// or DHP falls below the threshold are disabled.
const double test_threshold_rel = 0.05;
const double test_threshold_abs = 0.5;

const unsigned rnd_seed = 42;

string out_filename_base = "splicing_log";
string out_summary_base = "splicing_summary";

int iteration_width = max(static_cast<int>(log10(iterations)) + 1, 5);
int par_top_width_2 = max(iteration_width, 8);
const int par_top_width_1 = 26;
size_t name_width = 0;
const int p_width = 4 + 2;

// A Res object is created for each motif.
struct Res {
  vector<char> seqs[3];
  vector<array<unsigned, 3>> thread_perm_yes;
  vector<pair<double, double>> diff0_res[3];
  vector<double> diff1_res;
  streampos insert_pos0[3];
  streampos insert_pos1;
  pair<int, int> cnt_pos_neg0[3];
  pair<string, double> pred0[3];
  pair<int, int> cnt_pos_neg1;
  pair<string, double> pred1;
  bool disable0[3];
  bool disable1;
  static int disabled_total;

  Res() : disable0{false, false, false}, disable1(false) {
    for (int ii = 0; ii < 3; ++ii) {
      diff0_res[ii] = vector<pair<double, double>>(iterations + 1);
      seqs[ii].resize(no_seqs[ii]);
    }
    diff1_res = vector<double>(iterations + 1);
  }
};

int Res::disabled_total = 0;

map<string, Res> res;

struct {
  int count = no_threads;
  int chunk_size = -1;
  int chunk_remainder = -1;
  vector<int> seq_start;
  vector<thread> threads;
  vector<ofstream> outfiles;
  vector<string> filenames;
} threads_;

string get_search_params(int i, string logfile_base) {
  throw runtime_error("Not supported.");
}

int main(int argc, char** argv) {
  if (argc == 1 || argc == 2 && (0 == strcmp(argv[1], "-h") ||
                                 0 == strcmp(argv[1], "--help"))) {
    cout << "Usage: " << endl;
    cout << "  ./permute repeat_id [threads=" << threads_.count
         << "] [permutations=" << iterations << "]" << endl;
    return 0;
  }

  int out_id = stoi(argv[1]);
  
  if (argc > 2) threads_.count = stoi(argv[2]);
  if (argc > 3) {
    iterations = stoi(argv[3]);
    iteration_width = max(static_cast<int>(log10(iterations)) + 1, 5);
    par_top_width_2 = max(iteration_width, 8);
  }

  string logfile_base = "found_elements_";
  logfile_base = logfile_base + to_string(out_id) + "_";

  string run_time_strs[4];
  run_time_str(0);

  out_filename_base += "_" + to_string(out_id) + (repeat_run ? "a" : "");
  out_summary_base += "_" + to_string(out_id) + (repeat_run ? "a" : "");

  outfile.open(out_filename_base + ".txt");
  outfile << string(67, '#') << " Parameters" << endl << endl;
  outfile << left << setw(par_top_width_1) << "iterations " << right
          << setw(par_top_width_2) << iterations << endl;
  outfile << left << setw(par_top_width_1) << "p_threshold " << right
          << setw(par_top_width_2) << fixed << setprecision(2) << p_threshold
          << endl;
  outfile << left << setw(par_top_width_1) << "test_threshold_rel " << right
          << setw(par_top_width_2) << fixed << setprecision(2)
          << test_threshold_rel << endl;
  outfile << left << setw(par_top_width_1) << "test_threshold_abs " << right
          << setw(par_top_width_2) << fixed << setprecision(2)
          << test_threshold_abs << endl;
  outfile << left << setw(par_top_width_1) << "threads" << right
          << setw(par_top_width_2) << threads_.count << endl;
  outfile << left << setw(par_top_width_1 - 1) << "search_params "
          << (repeat_run ? "REPEAT_RUN" : get_search_params(0, logfile_base))
          << endl
          << endl;

  time_t now = chrono::system_clock::to_time_t(chrono::system_clock::now());
  tm tm_now{};
#ifdef _WIN32
  localtime_s(&tm_now, &now);
#else
  localtime_r(&now, &tm_now);
#endif
  outfile << left << setw(par_top_width_1) << "Time"
          << put_time(&tm_now, "%Y-%m-%dT%H:%M:%S") << endl;

  ostringstream run_times[3];

  pnt_top_summary({}, {}, {});

  outfile << string(67, '#') << " Summary" << endl << endl;

  for (int i = 0; i < 3; ++i) {
    ifstream data("motif_locations/" + logfile_base + group_str(i) + ".txt",
                  ios::binary);

    if (repeat_run) {
      data.seekg(-2, ios_base::end);

      for (;;) {
        char ch;
        if ((int)data.tellg() <= 0)
          throw runtime_error(
              "Error parsing found_elements from previous run.");
        data.get(ch);
        if (ch == '\n') break;
        data.seekg(-2, ios_base::cur);
      }

      string lastline;
      getline(data, lastline);
      no_seqs[i] = stoi(lastline) + 1;
      cout << "Parsed number of sequences in group " << group_str(i)
           << " from found_elements: " << no_seqs[i] << endl;
      data.seekg(0, ios::beg);
    }

    for (auto& kv : res) kv.second.seqs[i].resize(no_seqs[i]);

    string line;
    vector<vector<string>> tsv;
    for (; getline(data, line);) {
      stringstream lineStream(line);
      string cell;
      vector<string> parsedRow;
      for (; getline(lineStream, cell, '\t');) parsedRow.push_back(cell);

      tsv.push_back(parsedRow);
    }
    data.close();

    vector<vector<string>>::iterator row;
    vector<string>::iterator col;
    ptrdiff_t id_seq = 0, id_name = 0;

    for (row = tsv.begin(); row != tsv.end(); ++row) {
      col = row->begin();
      if (row == tsv.begin()) {
        for (; col != row->end(); ++col) {
          if (*col == "sequence_id") id_seq = col - row->begin();
          if (*col == "element_name") id_name = col - row->begin();
        }
        continue;
      }
      if (stoul(col[id_seq]) >= no_seqs[i])
        throw runtime_error("Error parsing number of sequences (" +
                            to_string(no_seqs[i]) + ").");

      if (res.find(col[id_name]) == res.end()) {
        res[col[id_name]] = Res{};
        if (col[id_name].length() > name_width)
          name_width = col[id_name].length();
      }
      res[col[id_name]].seqs[i][stoul(col[id_seq])] = true;
    }

    outfile << "BASE:" << group_str(i) << endl;

    pnt_res(i, 0, -1);
  }

  run_time_strs[1] = run_time_str(1);
  cout << "--- " << run_time_strs[1] << endl;

  outfile << endl << string(67, '#') << " Permutations" << endl << endl;

  threads_.chunk_size = iterations / threads_.count;
  threads_.chunk_remainder = iterations % threads_.count;

  threads_.seq_start.resize(threads_.count + 1);
  threads_.threads.reserve(threads_.count);
  threads_.outfiles.resize(threads_.count);
  threads_.filenames.resize(threads_.count);

  for (auto& kv : res) kv.second.thread_perm_yes.resize(threads_.count);

  for (int i = 0; i <= threads_.count; ++i)
    threads_.seq_start[i] =
        i * threads_.chunk_size + min(i, threads_.chunk_remainder);

  for (int i = 0; i < threads_.count; ++i) {
    threads_.filenames[i] = out_filename_base + "_" + to_string(i) + ".txt";
    threads_.outfiles[i].open(threads_.filenames[i]);
  }

  assert(threads_.seq_start[threads_.count] == iterations);

  mt19937 g(rnd_seed);
  uniform_int_distribution<unsigned> distr(0,
                                           (numeric_limits<unsigned>::max)());

  for (int i = 0; i < threads_.count; ++i)
    threads_.threads.emplace_back(it_step, threads_.seq_start[i],
                                  threads_.seq_start[i + 1], i, distr(g));

  for (int i = 0; i < threads_.count; ++i) threads_.threads[i].join();

  run_time_strs[2] = run_time_str(2);
  cout << "--- " << run_time_strs[2] << endl;

  for (int i = 0; i < threads_.count; ++i) {
    threads_.outfiles[i].close();
    ifstream infile(threads_.filenames[i]);
    outfile << infile.rdbuf();
    infile.close();
    remove(threads_.filenames[i].c_str());
  }

  int tests_cnt = static_cast<int>(res.size()) * 3 - Res::disabled_total;
  double adj_p = p_threshold / tests_cnt;

  for (int i = 2; i >= 0; --i) {
    for (auto rit = res.rbegin(); rit != res.rend(); ++rit) {
      auto& r = rit->second;

      if (i == 2 && !r.disable1) {
        auto& w = r.diff1_res;

        r.cnt_pos_neg1.first = static_cast<int>(
            count_if(w.begin() + 1, w.end(), [](auto dp) { return dp > 0; }));
        r.cnt_pos_neg1.second = static_cast<int>(
            count_if(w.begin() + 1, w.end(), [](auto dp) { return dp < 0; }));

        outfile.seekp(r.insert_pos1);
        outfile << "    <rnd " << right << setw(iteration_width)
                << r.cnt_pos_neg1.first << "    >rnd " << right
                << setw(iteration_width) << r.cnt_pos_neg1.second;

        if (w[0] < 0.0)
          r.pred1 = make_pair(
              ">", r.cnt_pos_neg1.second / static_cast<double>(iterations));
        else if (w[0] > 0.0)
          r.pred1 = make_pair(
              "<", r.cnt_pos_neg1.first / static_cast<double>(iterations));
        else
          r.pred1 = make_pair("?", 1.0);
      }

      if (i == 1 && !r.disable0[1] || i == 2 && !r.disable0[2]) {
        auto& v = r.diff0_res[i];

        r.cnt_pos_neg0[i].first = static_cast<int>(count_if(
            v.begin() + 1, v.end(), [](auto dp) { return dp.second > 0; }));
        r.cnt_pos_neg0[i].second = static_cast<int>(count_if(
            v.begin() + 1, v.end(), [](auto dp) { return dp.second < 0; }));

        outfile.seekp(r.insert_pos0[i]);
        outfile << "    <rnd " << right << setw(iteration_width)
                << r.cnt_pos_neg0[i].first << "    >rnd " << right
                << setw(iteration_width) << r.cnt_pos_neg0[i].second;

        if (v[0].second < 0.0)
          r.pred0[i] = make_pair(
              ">", r.cnt_pos_neg0[i].second / static_cast<double>(iterations));
        else if (v[0].second > 0.0)
          r.pred0[i] = make_pair(
              "<", r.cnt_pos_neg0[i].first / static_cast<double>(iterations));
        else
          r.pred0[i] = make_pair("?", 1.0);
      }

      if (i == 1 && !r.disable1) {
        outfile.seekp(r.insert_pos0[i]);
        outfile.seekp(2 * (strlen("    <rnd ") + iteration_width),
                      ios_base::cur);
        outfile << "    " << r.pred1.first << group_str(2) << " p=" << right
                << setw(p_width) << fixed << setprecision(p_width - 2)
                << r.pred1.second << (r.pred1.second <= adj_p ? '*' : ' ');
      }

      if (i == 0) {
        for (int ii = 1; ii < 3; ++ii)
          if (!r.disable0[ii]) {
            outfile.seekp(r.insert_pos0[i]);
            if (ii == 2)
              outfile.seekp(strlen("    >UHP p=*") + p_width, ios_base::cur);
            outfile << "    " << r.pred0[ii].first << group_str(ii)
                    << " p=" << right << setw(p_width) << fixed
                    << setprecision(p_width - 2) << r.pred0[ii].second
                    << (r.pred0[ii].second <= adj_p ? '*' : ' ');
          }
      }
    }
  }

  ofstream summary;
  summary.open(out_summary_base + ".txt");
  summary
      << "name\tt0.perc\tt1.perc\tt2.perc\tt0.vs.t1\tt0.vs.t2\tt1.vs.t2\tt0.vs."
         "t1.p\tt0.vs.t2.p\tt1.vs.t2.p\tt0.vs.t1.s\tt0.vs.t2.s\tt1.vs.t2.s"
      << endl;

  for (auto it = res.begin(); it != res.end(); ++it) {
    auto& r = it->second;
    summary << it->first << '\t';
    for (int i = 0; i < 3; ++i) summary << r.diff0_res[i][0].first << '\t';
    for (int i = 1; i < 3; ++i)
      summary << (r.disable0[i] ? "" : r.pred0[i].first)
              << (r.disable0[i] ? "" : group_str(i)) << '\t';
    summary << (r.disable1 ? "" : r.pred1.first)
            << (r.disable1 ? "" : group_str(2)) << '\t';
    for (int i = 1; i < 3; ++i)
      summary << (r.disable0[i] ? "" : to_string(r.pred0[i].second)) << '\t';
    summary << (r.disable1 ? "" : to_string(r.pred1.second)) << '\t';
    for (int i = 1; i < 3; ++i)
      summary << (r.disable0[i]
                      ? ""
                      : (r.pred0[i].second <= adj_p ? "TRUE" : "FALSE"))
              << '\t';
    summary << (r.disable1 ? "" : (r.pred1.second <= adj_p ? "TRUE" : "FALSE"))
            << endl;
  }

  run_time_strs[3] = run_time_str(3);
  cout << "--- " << run_time_strs[3] << endl;

  pnt_top_summary(tests_cnt, adj_p, run_time_strs);

  cout << "/// " << left << setw(par_top_width_1)
       << "Logfile: " << out_filename_base + ".txt" << endl;
  cout << "/// " << left << setw(par_top_width_1)
       << "Summary file: " << out_summary_base + ".txt" << endl;

  outfile.close();
  summary.close();
}

void it_step(int i_begin, int i_end, int thread_no, unsigned my_seed) {
  mt19937 g(my_seed);

  for (int k = i_begin; k < i_end; ++k) {
    vector<pair<int, unsigned>> v;
    v.reserve(no_seqs[0] + no_seqs[1] + no_seqs[2]);
    for (int i = 0; i < 3; ++i)
      for (unsigned j = 0; j < no_seqs[i]; ++j) v.emplace_back(i, j);

    shuffle(v.begin(), v.end(), g);

    for (int i = 0; i < 3; ++i) {
      for (auto& kv : res) kv.second.thread_perm_yes[thread_no][i] = 0;

      auto v_begin = v.begin();
      for (int ii = i - 1; ii >= 0; --ii) v_begin += no_seqs[ii];
      auto v_end = v_begin + no_seqs[i];

      for (auto& kv : res) {
        if (i == 0 && kv.second.disable0[1] && kv.second.disable0[2]) continue;
        if (i == 1 && kv.second.disable0[1] && kv.second.disable1) continue;
        if (i == 2 && kv.second.disable0[2] && kv.second.disable1) continue;

        for (auto v_it = v_begin; v_it < v_end; ++v_it)
          if (kv.second.seqs[v_it->first][v_it->second])
            ++kv.second.thread_perm_yes[thread_no][i];
      }

      threads_.outfiles[thread_no] << "rnd(" << k << "):" << group_str(i)
                                   << endl;
      pnt_res(i, k + 1, thread_no);

      if (i == 2 && v_end != v.end())
        throw runtime_error("Still elements left in vector.");
    }

    threads_.outfiles[thread_no] << endl;
  }
}

string group_str(int i) {
  if (i == 0) return "0";
  if (i == 1) return "UHP";
  if (i == 2) return "DHP";
  assert(false);
  return "";
}

bool fexists(string fname) {
#ifdef _WIN32
  DWORD fileAttributes = GetFileAttributesA(fname.c_str());
  if (fileAttributes == INVALID_FILE_ATTRIBUTES) return false;
  return !(fileAttributes & FILE_ATTRIBUTE_DIRECTORY);
#else
  return access(fname.c_str(), F_OK) != -1;
#endif
}

string run_time_str(int no) {
  ostringstream oss;
  double duration = 0;

#ifdef _WIN32
  static LARGE_INTEGER start, end, frequency;

  if (no > 0) {
    QueryPerformanceCounter(&end);
    duration =
        static_cast<double>(end.QuadPart - start.QuadPart) / frequency.QuadPart;
  }
#else
  static chrono::high_resolution_clock::time_point start, end;
  if (no > 0) {
    if (no >= 0) end = chrono::high_resolution_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
    duration = elapsed.count() / 1000.0;
  }
#endif

  if (abs(no) == 1)
    oss << left << setw(par_top_width_1) << "Read motif locations: " << right
        << setw(8) << fixed << setprecision(2) << duration << " seconds";
  if (abs(no) == 2)
    oss << left << setw(par_top_width_1) << "Run permutations: " << right
        << setw(8) << fixed << setprecision(2) << duration << " seconds";
  if (abs(no) == 3)
    oss << left << setw(par_top_width_1) << "Compute results: " << right
        << setw(8) << fixed << setprecision(2) << duration << " seconds";

#ifdef _WIN32
  if (no == 0) QueryPerformanceFrequency(&frequency);
  if (no >= 0) QueryPerformanceCounter(&start);
#else
  if (no >= 0) start = chrono::high_resolution_clock::now();
#endif

  return oss.str();
}

void pnt_top_summary(int cnt_tests, double adj_p, string run_time_strs[]) {
  static streampos top_pos = 0;
  if (top_pos == 0) {
    string possible_str = "( -";
    top_pos = outfile.tellp();

    outfile << run_time_str(-1) << endl
            << run_time_str(-2) << endl
            << run_time_str(-3) << endl
            << endl;
    outfile << left << setw(par_top_width_1) << "Total elements " << right
            << setw(par_top_width_2) << '-' << endl;
    outfile << left << setw(par_top_width_1) << "Completed tests " << right
            << setw(par_top_width_2) << '-' << "    " << right
            << setw(par_top_width_2) << possible_str << " possible, " << right
            << setw(strlen("100.00")) << '-' << "%)" << endl;
    outfile << left << setw(par_top_width_1) << "Adjusted p " << right
            << setw(par_top_width_2) << '-' << endl
            << endl;
  } else {
    string possible_str = "( " + to_string(res.size() * 3);
    outfile.seekp(top_pos);
    outfile << run_time_strs[1] << endl
            << run_time_strs[2] << endl
            << run_time_strs[3] << endl
            << endl;
    outfile << left << setw(par_top_width_1) << "Total elements " << right
            << setw(par_top_width_2) << res.size() << endl;
    outfile << left << setw(par_top_width_1) << "Completed tests " << right
            << setw(par_top_width_2) << cnt_tests << "    " << right
            << setw(par_top_width_2) << possible_str << " possible, " << right
            << setw(strlen("100.00")) << fixed << setprecision(2)
            << cnt_tests * 100.0 / (res.size() * 3) << "%)" << endl;
    outfile << left << setw(par_top_width_1) << "Adjusted p " << right
            << setw(par_top_width_2) << scientific << setprecision(2) << adj_p
            << endl
            << endl;
  }
}

void pnt_res(int i, int pos, int thread_no) {
  ofstream& outfile_here = pos == 0 ? outfile : threads_.outfiles[thread_no];
  for (auto& kv : res) {
    auto& r = kv.second;
    if (r.disable0[1] && r.disable0[2] && r.disable1) continue;

    size_t size = pos == 0 ? count(r.seqs[i].begin(), r.seqs[i].end(), true)
                           : r.thread_perm_yes[thread_no][i];
    auto& diff0 = r.diff0_res[i][pos];
    diff0.first = 100.0 * size / no_seqs[i];
    outfile_here << "  " << left << setw(name_width + 4) << kv.first << ": "
                 << right << setw(4) << fixed << setprecision(1) << diff0.first;

    if (i > 0) {
      diff0.second = diff0.first - r.diff0_res[0][pos].first;
      outfile_here << "    (" << right << setw(6) << fixed << setprecision(1)
                   << diff0.second << ")";

      if (pos > 0) {
        diff0.second = diff0.second - r.diff0_res[i][0].second;
        outfile_here << "    " << showpos << right << setw(6) << fixed
                     << setprecision(1) << diff0.second
                     << resetiosflags(ios::showpos);
      }
    }

    if (pos == 0) {
      r.insert_pos0[i] = outfile_here.tellp();

      if (i > 0) {
        outfile_here << "    <rnd " << right << setw(iteration_width) << '-'
                     << "    >rnd " << right << setw(iteration_width) << '-';

        if (abs(diff0.second) < test_threshold_abs ||
            abs(diff0.second) / r.diff0_res[0][0].first < test_threshold_rel) {
          r.disable0[i] = true;
          ++Res::disabled_total;
        }
      }

      if (i == 0)
        for (int ii = 1; ii < 3; ++ii)
          outfile_here << "    " << ' ' << group_str(ii) << " p=" << right
                       << setw(p_width) << '-' << ' ';
    }

    auto& diff1 = r.diff1_res[pos];

    if (i == 2) {
      diff1 = diff0.first - r.diff0_res[1][pos].first;
      outfile_here << "    (" << right << setw(6) << fixed << setprecision(1)
                   << diff1 << ")";

      if (pos > 0) {
        diff1 = diff1 - r.diff1_res[0];
        outfile_here << "    " << showpos << right << setw(p_width) << fixed
                     << setprecision(1) << diff1 << resetiosflags(ios::showpos);
      }
    }

    if (pos == 0) {
      if (i == 2) {
        r.insert_pos1 = outfile_here.tellp();
        outfile_here << "    <rnd " << right << setw(iteration_width) << '-'
                     << "    >rnd " << right << setw(iteration_width) << '-';
        if (abs(diff1) < test_threshold_abs ||
            abs(diff1) / r.diff0_res[1][0].first < test_threshold_rel) {
          r.disable1 = true;
          ++Res::disabled_total;
        }
      }

      if (i == 1)
        outfile_here << "    " << ' ' << group_str(2) << " p=" << right
                     << setw(p_width) << '-' << ' ';
    }

    outfile_here << endl;
  }
}
