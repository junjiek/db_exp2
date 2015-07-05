#include "SimJoiner.h"
#include <fstream>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <queue>
#include <cstring>

using namespace std;

inline int max3(int x, int y, int z) {
    return max(max(x, y), z);
}

inline int min3(int x, int y, int z) {
    return min(min(x, y), z);
}

// for debug
void print(vector<string>& v, int prefix_num) {
    for (int i = 0; i < prefix_num; i++)
        cout << v[i] << " ";
    cout << "| ";
    for (int i = prefix_num; i < (int)v.size(); i++)
        cout << v[i] << " ";
    cout << endl;
}

SimJoiner::SimJoiner() {}

SimJoiner::~SimJoiner() {}

void SimJoiner::readDataED(const char *filename1, const char *filename2) {
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

void SimJoiner::filterED(unsigned tau) {
    candidatesED.clear();
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
                                candidatesED.push_back(make_pair(id1, id2));
                                set.insert(id1);
                            }
                        }
                }
            }
        }
        for (auto& id1 : smallstr) {
            candidatesED.push_back(make_pair(id1, id2));
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

bool EDJoinResultCompare(const EDJoinResult &a, const EDJoinResult &b) {
    return a.id1 == b.id1 ? a.id2 < b.id2 : a.id1 < b.id1;
}

int SimJoiner::joinED(const char *filename1, const char *filename2,
                      unsigned threshold, vector<EDJoinResult> &result) {
    result.clear();
    readDataED(filename2, filename1);
    buildED(threshold);
    filterED(threshold);

    // verify
    vector<int> d0(maxLen1), d1(maxLen2);
    for(auto & i : candidatesED) {
        int ed = calDistED(str1[i.first], str2[i.second], threshold, d0, d1);
        if(ed <= (int)threshold)
            result.push_back(EDJoinResult(i.first, i.second, ed));
    }
    sort(result.begin(), result.end(), EDJoinResultCompare);

    return SUCCESS;
}


void SimJoiner::splitWords(char* str, vector<string>& words) {
    char* token = strtok(str, " ");
    unordered_set<string> exist;
    while (token != NULL) {
        string sub(token);
        if (exist.find(sub) == exist.end()) {
            exist.insert(sub);
            words.push_back(sub);
            wordListJac[sub] ++;
        }
        token = strtok(NULL, " ");
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
        for(int j = 0; j < prefix; ++ j) {
            string& sub = words1[id1][j];
            if(invertedListJac.find(sub) == invertedListJac.end()) {
                vector<int> v;
                v.push_back(id1);
                invertedListJac[sub] = v;
            } else {
                invertedListJac[sub].push_back(id1);
            }
        }
    }
}

int SimJoiner::calIntersec(int id1, int len1, int id2, int len2) {
    int i = 0, j = 0;
    int intersec = 0;

    while (i < len1 && j < len2) {
        if (words1[id1][i] == words2[id2][j]) {
            lastMatchPos.first = i;
            lastMatchPos.second = j;
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

void SimJoiner::filterJac(double tau, vector<JaccardJoinResult> &result) {
    candidatesJac.clear();
    for (int id2 = 0; id2 < (int)words2.size(); id2++) {
        int length = (int)words2[id2].size();
        int prefix = floor((1 - tau) * length + 1);
        if (length <= prefix) {
            for(int j = 0; j < (int)words1.size(); ++ j)
                candidatesJac.push_back(make_pair(j, id2));
            continue;
        }
        unordered_set<int> candidates;
        for(int j = 0; j < prefix; ++ j) {
            string& sub = words2[id2][j];
            for(auto & id1 : invertedListJac[sub]) {
                if (candidates.find(id1) != candidates.end())
                    continue;
                else
                    candidates.insert(id1);
                int length1 = words1[id1].size();
                if (length1 >= tau * length && length >= tau * length1) {
                    
                    int threshold = ceil(tau * (length + length1) / (1 + tau));
                    int prefix1 = length1 - threshold + 1;
                    int prefix2 = length - threshold + 1;
                    int prefixIntersec = calIntersec(id1, prefix1, id2, prefix2);
                    int ubound = prefixIntersec + min(length - prefix2, length1 - prefix1);
                    if (ubound >= threshold) {
                        prefix1 = min(prefix1, 1 + lastMatchPos.first);
                        prefix2 = min(prefix2, 1 + lastMatchPos.second);
                        int suffixIntersec =
                            calSuffixIntersec(words1[id1], make_pair(prefix1, length1),
                                              words2[id2], make_pair(prefix2, length),
                                              threshold - prefixIntersec);
                        int intersec = prefixIntersec + suffixIntersec;
                        if (intersec >= threshold) {
                            double jac = (double) intersec / (words1[id1].size() + words2[id2].size() - intersec);
                            result.push_back(JaccardJoinResult(id1, id2, jac));
                        }
                    }
                }
            }
        }
        for (auto& id1 : smallstr) {
            candidatesJac.push_back(make_pair(id1, id2));
        }
    }

}

// int SimJoiner::calUbound(vector<string>& words1, pair<int, int> pos1,
//                          vector<string>& words2, pair<int, int> pos2, int depth) {
//     int len1 = pos1.second - pos1.first;
//     int len2 = pos2.second - pos2.first;
//     if (depth == 0 || len1 == 0 || len2 == 0) {
//         return min(len1, len2);
//     }
//     int mid1 = (pos1.first + pos1.second) >> 1;
//     int mid2 = binSearch(words2, words1[mid1], pos2.first, pos2.second);
//     // cout << words1[mid1] << mid1 << ", " << words2[mid2] << mid2 << endl;
//     return calUbound(words1, make_pair(pos1.first, mid1 + 1),
//                      words2, make_pair(pos2.first, mid2 + 1), depth - 1)
//            + calUbound(words1, make_pair(mid1 + 1, pos1.second),
//                        words2, make_pair(mid2 + 1, pos2.second), depth - 1);
// }

int intervalLen(const interval & i) {
    return i.second - i.first;
}
int intervalPairLen(const intervalPair & ip) {
    return min(intervalLen(ip.first), intervalLen(ip.second));
}

struct intervalPairLenCmp {
    bool operator() (const intervalPair& a, const intervalPair& b) {
        return intervalPairLen(a) < intervalPairLen(b);
    }
};


int SimJoiner::calSuffixIntersec(vector<string>& words1, pair<int, int> pos1,
                                 vector<string>& words2, pair<int, int> pos2,
                                 int threshold) {
    priority_queue<intervalPair, vector<intervalPair>, intervalPairLenCmp> pq;
    intervalPair ip1 = make_pair(pos1, pos2);
    int ubound = intervalPairLen(ip1);
    pq.push(ip1);
    int intersec = 0;
    while (!pq.empty()) {
        intervalPair ip = pq.top();
        pq.pop();
        ubound -= intervalPairLen(ip);
        int l1 = (ip.first).first;
        int r1 = (ip.first).second;
        int mid1 = (l1 + r1) >> 1;
        int l2 = (ip.second).first;
        int r2 = (ip.second).second;
        int mid2 = binSearch(words2, words1[mid1], l2, r2);
        bool match = words1[mid1] == words2[mid2];
        if (match) {
            intersec ++;
            ubound ++;
        }
        intervalPair left = make_pair(make_pair(l1, mid1), make_pair(l2, match ? mid2 : mid2+1));
        intervalPair right = make_pair(make_pair(mid1+1, r1), make_pair(mid2+1, r2));
        int leftLen = intervalPairLen(left);
        int rightLen = intervalPairLen(right);
        if (leftLen > 0) {
            pq.push(left);
            ubound += leftLen;
        }
        if (rightLen > 0) {
            pq.push(right);
            ubound += rightLen;
        }
        if (ubound < threshold)
            return 0;
    }
    return intersec;
}

int SimJoiner::binSearch(vector<string> & words, string& sub, int lo, int hi) {
    while (lo  < hi) {
        int mi = (lo + hi) >> 1;
        string& mid_sub = words[mi];
        if (wordListJac[sub] < wordListJac[mid_sub] || 
            (wordListJac[sub] == wordListJac[mid_sub] && sub < mid_sub))
            hi = mi;
        else
            lo = mi + 1;
    }
    return --lo;
}

void SimJoiner::readDataJac(const char *filename1, const char *filename2) {
    size_t len = 0;
    char *line = NULL;
    words1.clear();
    words2.clear();
    FILE *fp = fopen(filename1, "r");
    while (getline(&line, &len, fp) != -1) {
        vector<string> wordList;
        line[strlen(line)-1] = '\0';
        splitWords(line, wordList);
        words1.push_back(wordList);
    }
    fclose(fp);
    fp = fopen(filename2, "r");
    while (getline(&line, &len, fp) != -1) {
        vector<string> wordList;
        line[strlen(line)-1] = '\0';
        splitWords(line, wordList);
        words2.push_back(wordList);
    }
}

// sort by frequency
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

bool JaccardJoinResultCompare(const JaccardJoinResult &a, const JaccardJoinResult &b) {
    return a.id1 == b.id1 ? a.id2 < b.id2 : a.id1 < b.id1;
}

double SimJoiner::calDistJac(int id1, int id2, double threshold) {
    int dataSize = (int)words1[id1].size();
    int querySize = (int)words2[id2].size();
    if (dataSize * threshold > querySize
        || querySize * threshold > dataSize)
        return 0;
    int intersec = calIntersec(id1, dataSize, id2, querySize);
    return (double)intersec / (dataSize + querySize - intersec);
}

int SimJoiner::joinJaccard(const char *filename1, const char *filename2,
                           double threshold, vector<JaccardJoinResult> &result) {
    result.clear();
    // start = clock();
    readDataJac(filename2, filename1);
    // finish = clock();
    // printf("readDataJac: %lfs \n", (double)(finish - start) / CLOCKS_PER_SEC);

    start = clock();    
    for (auto & v : words1) {
        qsort(v, 0, v.size());
    }
    for (auto & v : words2) {
        qsort(v, 0, v.size());
    }
    // finish = clock();
    // printf("sort: %lfs \n", (double)(finish - start) / CLOCKS_PER_SEC);
    
    // start = clock();
    buildJac(threshold);
    // finish = clock();
    // printf("buildJac: %lfs \n", (double)(finish - start) / CLOCKS_PER_SEC);

    // start = clock();
    filterJac(threshold, result);
    // finish = clock();
    // printf("filterJac: %lfs \n", (double)(finish - start) / CLOCKS_PER_SEC);
    
    // verify
    // start = clock();
    // cout << candidatesJac.size() << endl;
    for(auto & i : candidatesJac) {
        double jac = calDistJac(i.first, i.second, threshold);
        if(jac >= threshold - EPS)
            result.push_back(JaccardJoinResult(i.first, i.second, jac));
    }
    // finish = clock();
    // printf("verify: %lfs \n", (double)(finish - start) / CLOCKS_PER_SEC);
    
    sort(result.begin(), result.end(), JaccardJoinResultCompare);

    return SUCCESS;
}

