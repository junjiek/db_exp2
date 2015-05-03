#include "SimJoiner.h"
#include <fstream>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <unordered_set>
using namespace std;

bool EDJoinResultCompare(const EDJoinResult &a, const EDJoinResult &b) {
    return a.id1 == b.id1 ? a.id2 < b.id2 : a.id1 < b.id1;
}

bool JaccardJoinResultCompare(const JaccardJoinResult &a, const JaccardJoinResult &b) {
    return a.id1 == b.id1 ? a.id2 < b.id2 : a.id1 < b.id1;
}

SimJoiner::SimJoiner() {
    q = 4;
}

SimJoiner::~SimJoiner() {}

void SimJoiner::buildJac(double threshold) {
    smallstr.clear();
    invertedListJac.clear();
    for (int i = 0; i < (int)str1.size(); i++) {
        int length = str1[i].length();
        int prefix = floor((1 - threshold) * length + 1);
        if (length < prefix + (int)q - 1) {
            smallstr.push_back(i);
            continue;
        }
        for(int j = 0; j < prefix; ++ j) {
            string sub = str1[i].substr(j, q);
            if(invertedListJac.find(sub) == invertedListJac.end()) {
                vector<pair<unsigned, unsigned>> v;
                v.push_back(make_pair(i, j));
                invertedListJac[sub] = v;
            }
            else {
                invertedListJac[sub].push_back(make_pair(i, j));
            }
        }
    }
}


void SimJoiner::filterJac(double threshold) {
    for (int id2 = 0; id2 < (int)str2.size(); id2++) {
        int length = str2[id2].length();
        int prefix = length - floor(threshold * length) + 1;
        if(length < prefix + (int)q - 1) {
            for(int j = 0; j < (int)str1.size(); ++ j)
                rawResultJac.push_back(JaccardJoinResult(j, id2, 0));
            continue;
        }
        unordered_map<int, int> m;
        for(int j = 0; j < prefix; ++ j) {
            string sub = str2[id2].substr(j, q);
            for(auto & p : invertedListJac[sub]) {
                int length2 = str1[p.first].length();
                if(length2 >= threshold * length) {
                    int alpha = floor(threshold * (length + length2) / (1 + threshold));
                    int ubound = min(length - j, length2 - int(p.second));
                    if(m[p.first] + ubound >= alpha)
                        m[p.first] += 1;
                }
            }
        }
        for (auto & p : m) {
            if(p.second > 0)
                rawResultJac.push_back(JaccardJoinResult(p.first, id2, p.second));
        }
        for (auto& ss : smallstr) {
            rawResultJac.push_back(JaccardJoinResult(ss, id2, 0));
        }
    }
}

double SimJoiner::calDistJac(string &a, string &b, int T) {
    double dis = 0;
    unsigned len_a = a.length(), len_b = b.length();
    if(len_a < q || len_b < q)
        return dis;
    unordered_map<string, int> m1, m2;
    for(unsigned i = 0; i <= len_a - q; ++ i) {
        string sub = a.substr(i, q);
        m1[sub] += 1;
    }
    for(unsigned i = 0; i <= len_b - q; ++ i) {
        string sub = b.substr(i, q);
        m2[sub] += 1;
    }
    T = 0;
    for(auto & p : m1)
        T += min(p.second, m2[p.first]);
    dis = (double)T / (len_a + len_b - 2 * (q - 1) - T);
    return dis;
}

int SimJoiner::joinJaccard(const char *filename1, const char *filename2,
                           double threshold, vector<JaccardJoinResult> &result) {
    result.clear();
    readData(filename2, filename1);
    buildJac(threshold);
    filterJac(threshold);

    // verify
    for(auto & i : rawResultJac) {
        double jac = calDistJac(str1[i.id1], str2[i.id2], i.s);
        if(jac >= threshold)
            result.push_back(JaccardJoinResult(i.id1, i.id2, jac));
    }
    sort(result.begin(), result.end(), JaccardJoinResultCompare);

    return SUCCESS;
}

void SimJoiner::insertED(unsigned id, unsigned length, unsigned pos, string &str) {
    if(invertedListED.size() < length + 1)
        invertedListED.resize(length + 1);
    if(invertedListED[length].size() < pos + 1)
        invertedListED[length].resize(pos + 1);
    if(invertedListED[length][pos].find(str) == invertedListED[length][pos].end()) {
        vector<int> v;
        v.push_back(id);
        invertedListED[length][pos][str] = v;
    } else {
        invertedListED[length][pos][str].push_back(id);
    }

}

void SimJoiner::buildED(unsigned tau) {
    smallstr.clear();
    for(int i = 0; i < (int)str1.size(); ++ i) {
        unsigned length = str1[i].length();
        if(length < tau + 1) {
            smallstr.push_back(i);
            continue;
        }
        // divide into tau+1 groups
        int shortLen = length / (tau + 1), longLen = shortLen + 1;
        int longNum = length % (tau + 1), shortNum = (int)tau + 1 - longNum; 
        for(int j = 0; j < shortNum; ++ j) {
            string sub = str1[i].substr(j * shortLen, shortLen);
            insertED(i, length, j, sub);
        }
        for(int j = 0; j < longNum; ++ j) {
            string sub = str1[i].substr(shortNum * shortLen + j * longLen, longLen);
            insertED(i, length, shortNum + j, sub);
        }
    }
}

inline int max3(int x, int y, int z) {
    return max(max(x, y), z);
}

inline int min3(int x, int y, int z) {
    return min(min(x, y), z);
}

void SimJoiner::filterED(unsigned tau) {
    for(int id2 = 0; id2 < (int)str2.size(); ++ id2) {
        unsigned length = str2[id2].length();
        unordered_set<int> set;
        int minLen = max(1, int(length) - int(tau));
        int maxLen = min((int)(length + tau), (int)invertedListED.size() - 1);
        for (int len = minLen; len <= maxLen; ++ len) {
            for (int j = 0; j < (int)invertedListED[len].size(); j++) {
                if (invertedListED[len][j].size() == 0)
                    continue;
                int pos, substrlen;
                int longNum = len % (tau + 1), shortNum = (int)tau + 1 - longNum;
                int shortLen = len / (tau + 1), longLen = shortLen + 1;
                if (j >= shortNum) {
                    substrlen = longLen;
                    pos = shortLen * shortNum +  longLen * (j - shortNum) ;
                } else {
                    substrlen = shortLen;
                    pos = shortLen * j;
                }
                int delta = int(length) - len;
                // Multi-match-aware method
                int start = max3(0, pos - j, pos + delta - (int)tau + j - 2);
                int end = min3(int(length) - substrlen, pos + j, pos + delta + (int)tau - j);
                for(int k = start; k <= end; ++ k) {
                    string sub = str2[id2].substr(k, substrlen);
                    if(invertedListED[len][j].find(sub) != invertedListED[len][j].end())
                        for(auto & id1 : invertedListED[len][j][sub]) {
                            if(set.find(id1) == set.end()) {
                                rawResultED.push_back(make_pair(id1, id2));
                                set.insert(id1);
                            }
                        }
                }
            }
        }
        for (auto& id1 : smallstr) {
            rawResultED.push_back(make_pair(id1, id2));
        }
    }
}

unsigned SimJoiner::calDistED(string &a, string &b, unsigned threshold,
                              vector<int> &d0, vector<int> &d1) {
    unsigned dis = threshold + 1;
    int len_a = a.length(), len_b = b.length();
    if(abs(len_a - len_b) > threshold)
        return dis;
    for(int i = 0; i <= len_a; ++ i) {
        int l = max(0, i - (int)threshold), r = min(len_b, i + (int)threshold);
        int minDis = threshold + 1;
        for(int j = l; j <= r; ++ j) {
            if(i == 0)
                d1[j] = j;
            else if(j == 0)
                d1[j] = i;
            else {
                if(a[i - 1] == b[j - 1])
                    d1[j] = d0[j - 1];
                else
                    d1[j] = d0[j - 1] + 1;
                if(j > l) d1[j] = min(d1[j], d1[j - 1] + 1);
                if(j < i + (int)threshold) d1[j] = min(d1[j], d0[j] + 1);
            }
            minDis = min(minDis, d1[j]);    
        }
        if(minDis > (int)threshold)
            return dis;
        swap(d0, d1);
    }
    dis = d0[len_b];
    return dis;
}


int SimJoiner::joinED(const char *filename1, const char *filename2,
                      unsigned threshold, vector<EDJoinResult> &result) {
    result.clear();
    readData(filename2, filename1);
    buildED(threshold);
    filterED(threshold);

    // verify
    vector<int> d0(maxLen1), d1(maxLen2);
    for(auto & i : rawResultED) {
        int ed = calDistED(str1[i.first], str2[i.second], threshold, d0, d1);
        if(ed <= (int)threshold)
            result.push_back(EDJoinResult(i.first, i.second, ed));
    }
    sort(result.begin(), result.end(), EDJoinResultCompare);

    return SUCCESS;
}

void SimJoiner::readData(const char *filename1, const char *filename2) {
    ifstream fin1(filename1);
    string line;
    maxLen1 = 0;
    while(getline(fin1, line)) {
        str1.push_back(line);
        if(line.length() > maxLen1)
            maxLen1 = line.length();
    }
    fin1.close();

    ifstream fin2(filename2);
    maxLen2 = 0;
    while(getline(fin2, line)) {
        str2.push_back(line);
        if(line.length() > maxLen2)
            maxLen2 = line.length();
    }
    fin2.close();
}
