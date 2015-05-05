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
    q = 10;
}

SimJoiner::~SimJoiner() {}


void SimJoiner::splitWords(string & str, vector<string>& words, bool buildWordList) {
    int nend = 0; 
    int nbegin = 0;
    unordered_set<string> exist;
    string sub;
    while (nend != -1) {
        nend = str.find_first_of(' ', nbegin);
        if (nend == -1) {
            sub = str.substr(nbegin, str.length() - nbegin);
        } else {
            sub = str.substr(nbegin, nend - nbegin);
        }
        if (exist.find(sub) == exist.end()) {
            exist.insert(sub);
            words.push_back(sub);
            if (buildWordList)
                wordListJac[sub] ++;
        }
        nbegin = nend + 1;
    }
}

void SimJoiner::buildJac(double threshold) {
    smallstr.clear();
    invertedListJac.clear();
    for (int id1 = 0; id1 < (int)words1.size(); id1++) {
        int length = (int)words1[id1].size();
        int prefix = floor((1 - threshold) * length + 1);
        if (length <= prefix) {
            smallstr.push_back(id1);
            continue;
        }
        // cout << id1 << " " << prefix << ": ";
        for(int j = 0; j < prefix; ++ j) {
            string& sub = words1[id1][j];
            // cout << sub << " ";
            if(invertedListJac.find(sub) == invertedListJac.end()) {
                vector<pair<unsigned, unsigned>> v;
                v.push_back(make_pair(id1, j));
                invertedListJac[sub] = v;
            } else {
                invertedListJac[sub].push_back(make_pair(id1, j));
            }
        }
        // cout << "| ";
        // for (int j = prefix; j < length; j++)
        //     cout << words1[id1][j] << " ";
        // cout << endl;
    }
    // for (auto & list : invertedListJac) {
    //     cout << list.first << ": ";
    //     for (auto & p : list.second) {
    //         cout << "(" <<  p.first << ", " << p.second << ") ";
    //     }
    //     cout << endl;
    // }
}

int SimJoiner::calIntersec(int id1, int len1, int id2, int len2) {
    int i = 0, j = 0;
    int intersec = 0;
    while (i < len1 && j < len2) {
        if (words1[id1][i] == words2[id2][j]) {
            intersec ++;
            i ++;
            j ++;
        } else if (wordListJac[words1[id1][i]] > wordListJac[words2[id2][j]] || 
                   (wordListJac[words1[id1][i]] == wordListJac[words2[id2][j]] && 
                    words1[id1][i] > words2[id2][j])) {
            j ++;
        } else {
            i ++;
        }
    }
    return intersec;
}

void SimJoiner::filterJac(double tau) {
    rawResultJac.clear();
    for (int id2 = 0; id2 < (int)words2.size(); id2++) {
        int length = (int)words2[id2].size();
        int prefix = floor((1 - tau) * length + 1);
        if (length <= prefix) {
            for(int j = 0; j < (int)words1.size(); ++ j)
                rawResultJac.push_back(make_pair(j, id2));
            continue;
        }
        unordered_map<int, int> m;
        for(int j = 0; j < prefix; ++ j) {
            string& sub = words2[id2][j];
            for(auto & p : invertedListJac[sub]) {
                int length2 = words1[p.first].size();
                if (length2 >= tau * length) {
                    int threshold = ceil(tau * (length + length2) / (1 + tau));
                    int ubound = min(length - j, length2 - int(p.second));
                    if (m[p.first] + ubound >= threshold)
                        m[p.first] += 1;

                }
            }
        }
        for (auto & p : m) {
            if(p.second > 0)
                rawResultJac.push_back(make_pair(p.first, id2));
        }
        // unordered_set<int> candidates;
        // for(int j = 0; j < prefix; ++ j) {
        //     string& sub = words2[id2][j];
        //     for(auto & id1 : invertedListJac[sub]) {
        //         if (candidates.find(id1.first) != candidates.end())
        //             continue;
        //         else
        //             candidates.insert(id1.first);
        //         int length1 = words1[id1.first].size();
        //         if (length1 >= tau * length) {
                    
        //             int prefix1 = floor((1 - tau) * length1 + 1);
        //             int prefixIntersec = calIntersec(id1.first, prefix1, id2, prefix);
        //             int ubound = prefixIntersec + min(length - prefix, length1 - prefix1);
        //             int threshold = ceil(tau * (length + length1) / (1 + tau));
        //             if (ubound >= threshold)
        //                 rawResultJac.push_back(make_pair(id1.first, id2));
        //         }
        //     }
        // }
        for (auto& ss : smallstr) {
            rawResultJac.push_back(make_pair(ss, id2));
        }
    }
}

double SimJoiner::calDistJac(int id1, int id2) {
    int intersec = calIntersec(id1, words1[id1].size(), id2, words2[id2].size());
    return (double)intersec / (words1[id1].size() + words2[id2].size() - intersec);
}


void SimJoiner::readDataJac(const char *filename1, const char *filename2) {
    ifstream fin1(filename1);
    string line;
    words1.clear();
    words2.clear();
    while(getline(fin1, line)) {
        vector<string> wordList;
        // cout << line << endl;
        splitWords(line, wordList, true);
        words1.push_back(wordList);
    }
    fin1.close();

    ifstream fin2(filename2);
    while(getline(fin2, line)) {
        vector<string> wordList;
        splitWords(line, wordList, false);
        words2.push_back(wordList);
    }
    fin2.close();
}


void SimJoiner::qsort(vector<string> & words, int lo, int hi) {
    if (hi - lo < 2) return;
    int mi = partition(words, lo, hi - 1);
    qsort(words, lo, mi);
    qsort(words, mi + 1, hi);
}

int SimJoiner::partition(vector<string> & words, int lo, int hi) {
    swap(words[lo], words[lo + rand() % (hi - lo + 1)]);
    string& pivot = words[lo];
    int mi = lo;
    for (int k = lo + 1; k <= hi; k ++) {
        if (wordListJac[words[k]] < wordListJac[pivot] || 
            (wordListJac[words[k]] == wordListJac[pivot] && words[k] < pivot) ) {
            swap(words[++mi], words[k]);
        }
    }
    swap(words[lo], words[mi]);
    return mi;
}

int SimJoiner::joinJaccard(const char *filename1, const char *filename2,
                           double threshold, vector<JaccardJoinResult> &result) {
    result.clear();
    start = clock();
    readDataJac(filename2, filename1);
    finish = clock();
    printf("readDataJac: %lfs \n", (double)(finish - start) / CLOCKS_PER_SEC);

    start = clock();    
    for (auto & v : words1) {
        qsort(v, 0, v.size());
    }
    for (auto & v : words2) {
        qsort(v, 0, v.size());
    }
    finish = clock();
    printf("sort: %lfs \n", (double)(finish - start) / CLOCKS_PER_SEC);
    
    start = clock();
    buildJac(threshold);
    finish = clock();
    printf("buildJac: %lfs \n", (double)(finish - start) / CLOCKS_PER_SEC);

    start = clock();
    filterJac(threshold);
    finish = clock();
    printf("filterJac: %lfs \n", (double)(finish - start) / CLOCKS_PER_SEC);
    
    // verify
    start = clock();
    cout << rawResultJac.size() << endl;
    for(auto & i : rawResultJac) {
        double jac = calDistJac(i.first, i.second);
        if(jac >= threshold)
            result.push_back(JaccardJoinResult(i.first, i.second, jac));
    }
    finish = clock();
    printf("verify: %lfs \n", (double)(finish - start) / CLOCKS_PER_SEC);
    
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
    invertedListED.clear();
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
    rawResultED.clear();
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
    str1.clear();
    str2.clear();
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
