#include "SimJoiner.h"
#include <cassert>
#include <unistd.h>
#include <err.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>

int main(int argc, char **argv)
{
  if (argc != 6) {
    fprintf(stderr, "Usage: %s [ed|j] threshold file1 file2 [0|1]\n", argv[0]);
    return 1;
  }
  size_t len = 0;
  char *line = NULL;
  vector<string> doc1, doc2;

  FILE *fp = fopen(argv[3], "r");
  if (! fp)
    errx(2, "%s not readable", argv[3]);
  while (getline(&line, &len, fp) != -1) {
    line[strlen(line)-1] = '\0';
    doc1.push_back(line);
  }

  fp = fopen(argv[4], "r");
  if (! fp)
    errx(2, "%s not readable", argv[4]);
  while (getline(&line, &len, fp) != -1) {
    line[strlen(line)-1] = '\0';
    doc2.push_back(line);
  }

  bool debug = atoi(argv[5]);

  SimJoiner joiner;
  string op = argv[1];
  clock_t start, finish;
  if (op == "j") {
    vector<JaccardJoinResult> res;
    float th = atof(argv[2]);
    start = clock();
    joiner.joinJaccard(argv[3], argv[4], th, res);
    finish = clock();
    printf("joinJaccard: %lfs passed, %d results found\n", (double)(finish - start) / CLOCKS_PER_SEC, (int)res.size());
    getchar();
    if (debug)
      for (auto i: res)
        printf("+ %lf\n  %s\n  %s\n", i.s, doc1[i.id1].c_str(), doc2[i.id2].c_str());
    else
      for (auto i: res)
        printf("+ %u %u %lf\n", i.id1, i.id2, i.s);
  } else if (op == "ed") {
    vector<EDJoinResult> res;
    unsigned th = atoi(argv[2]);
    start = clock();
    joiner.joinED(argv[3], argv[4], th, res);
    finish = clock();
    printf("joinEd: %lfs passed, %d results found\n", (double)(finish - start) / CLOCKS_PER_SEC, (int)res.size());
    getchar();
    if (debug)
      for (auto i: res)
        printf("+ %u\n  %s\n  %s\n", i.s, doc1[i.id1].c_str(), doc2[i.id2].c_str());
    else
      for (auto i: res)
        printf("+ %u %u %u\n", i.id1, i.id2, i.s);
  } else
    errx(3, "unknown operation");

  free(line);

  return 0;
}
