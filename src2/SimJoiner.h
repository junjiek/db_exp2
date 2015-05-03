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

const int SUCCESS = 0;
const int FAILURE = 1;

class SimJoiner {
private:
    unsigned maxLen1, maxLen2;
    // data from file1 and file2
    vector<string> str1, str2;
    vector<vector<string>> words1, words2;
    // invertedListED[length][pos][str] -> vector of id
    vector<vector<unordered_map<string, vector<int>>>> invertedListED;
    // str -> vector of pair<id, pos>
    unordered_map<string, vector<pair<unsigned, unsigned>>> invertedListJac;
    vector<unsigned> smallstr;
    vector<pair<unsigned, unsigned>> rawResultJac;
    vector<pair<unsigned, unsigned>> rawResultED;
    unsigned q;
    void insertED(unsigned id, unsigned length, unsigned pos, string &str);
    void readData(const char *filename1, const char *filename2);
    void readDataJac(const char *filename1, const char *filename2);
    void splitWords(string & str, vector<string>& words);
    void buildJac(double tau);
    void filterJac(double tau);
    double calDistJac(int id1, int id2);
    void buildED(unsigned threshold);
    void filterED(unsigned tau);
    unsigned calDistED(string &a, string &b, unsigned threshold,
                       vector<int> &d0, vector<int> &d1);
public:
    SimJoiner();
    ~SimJoiner();
    int joinJaccard(const char *filename1, const char *filename2,
                    double threshold, std::vector<JaccardJoinResult> &result);
    int joinED(const char *filename1, const char *filename2, unsigned threshold,
               std::vector<EDJoinResult> &result);
};
#endif
