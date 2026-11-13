// test_api.cpp - proves the C API works exactly as Kotlin/JNI will call it (text in, text out, no files).
// Build: g++ -std=c++20 -O2 -pthread -DEBOV_NO_MAIN test_api.cpp code.cpp -o test_api
extern "C" { const char* ebov_version(); char* ebov_run(const char*, int); char* ebov_screen(const char*, const char*, const char*); void ebov_free(char*); }
#include <cstdio>
#include <cstring>
int main() {
  std::printf("%s\n", ebov_version());
  const char* params = "tip_dose = 1e9\ntip_time = 1\nc_T = 0.6\nalpha = 0.05\nbeta_T = 2e-6\nV_ext = 2e-4\next_floor = 1\n";
  char* s = ebov_run(params, 0);   std::printf("--- summary ---\n%s", s); ebov_free(s);
  char* c = ebov_run(params, 1);   std::printf("--- csv: %zu bytes, first line: %.40s\n", std::strlen(c), c); ebov_free(c);
  char* v = ebov_run(params, 2);   std::printf("--- svg: %zu bytes, starts: %.30s\n", std::strlen(v), v); ebov_free(v);
  char* bad = ebov_run("alpha = 7\nnot_a_param = 3\n", 0); std::printf("--- bad input ---\n%s", bad); ebov_free(bad);
  char* r = ebov_screen("ext_floor = 1\n", "samples = 100\nc_T = 0.5 15 log\ntip_dose = 1e5 1e9 log\ntip_time = 0 5\nbeta_T = 2.5e-7 2e-6 log\n", "/tmp/api_screen");
  std::printf("--- screen report (first 3 lines) ---\n"); int n = 0; for (char* p = r; *p && n < 3; ++p) { std::putchar(*p); if (*p == '\n') ++n; } ebov_free(r);
  return 0;
}
