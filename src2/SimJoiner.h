#ifndef __EXP2_SIMJOINER_H__
#define __EXP2_SIMJOINER_H__

#include <vector>
#include <unordered_map>

using namespace std;

template <typename IDType, typename SimType>
struct JoinResult {
    JoinResult(IDType _id1, IDType _id2, SimType _s) : id1(_id1), id2(_id2), s(_s) {}
    IDType id1;
    IDType id2;
    SimType s;
};

typedef JoinResult<unsigned, double> JaccardJoinResult;
typedef JoinResult<unsigned, unsigned> EDJoinResult;
typedef pair<int, int> interval;
typedef pair<interval, interval> intervalPair;

const double EPS = 1e-8;
const int SUCCESS = 0;
const int FAILURE = 1;


class SimJoiner {
private:
    clock_t start, finish;
    unsigned maxLen1, maxLen2;
    pair<int, int> lastMatchPos;
    // data from file1 and file2
    vector<string> str1, str2;
    vector<vector<string>> words1, words2;
    // invertedListED[length][pos][str] -> vector of id
    vector<vector<unordered_map<string, vector<int>>>> invertedListED;
    // str -> vector of pair<id, pos>
    unordered_map<string, vector<int>> invertedListJac;
    unordered_map<string, int> wordListJac;
    vector<unsigned> smallstr;
    vector<pair<unsigned, unsigned>> candidatesJac;
    vector<pair<unsigned, unsigned>> candidatesED;
    void insertED(unsigned id, unsigned length, unsigned pos, string &str);
    void readDataED(const char *filename1, const char *filename2);
    void readDataJac(const char *filename1, const char *filename2);
    void splitWords(char* str, vector<string>& words);
    void buildJac(double tau);
    void filterJac(double tau, vector<JaccardJoinResult> &result);
    double calDistJac(int id1, int id2, double threshold);
    void buildED(unsigned threshold);
    void filterED(unsigned tau);
    unsigned calDistED(string &a, string &b, unsigned threshold,
                       vector<int> &d0, vector<int> &d1);
    int binSearch(vector<string> & words, string& sub, int lo, int hi);
    // int calUbound(vector<string>& words1, pair<int, int> pos1,
    //                      vector<string>& words2, pair<int, int> pos2, int depth);
    int calSuffixIntersec(vector<string>& words1, pair<int, int> pos1,
                          vector<string>& words2, pair<int, int> pos2,
                          int threshold);
    int partition(vector<string> & words, int lo, int hi);
    void qsort(vector<string> & words, int lo, int hi);
    int calIntersec(int id1, int len1, int id2, int len2);
public:
    SimJoiner();
    ~SimJoiner();
    int joinJaccard(const char *filename1, const char *filename2,
                    double threshold, std::vector<JaccardJoinResult> &result);
    int joinED(const char *filename1, const char *filename2, unsigned threshold,
               std::vector<EDJoinResult> &result);
};
#endif
