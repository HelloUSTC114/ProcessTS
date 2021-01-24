#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TChain.h"
#include "TSystem.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

enum DataType
{
    mppc,
    t0,
    undefined,
};
const int gLoopSize = 200;

struct LoopIndex
{
    UInt_t fIndex = 0;
    UInt_t fLoopSize = gLoopSize;

    LoopIndex(UInt_t size = gLoopSize) : fLoopSize(size) { fIndex = 0; }
    LoopIndex(UInt_t size, UInt_t index) { fLoopSize = (size), fIndex = (index % fLoopSize); }

    LoopIndex
    operator+(Int_t t)
    {
        LoopIndex a(fLoopSize, fIndex + t);
        return a;
    }
    // operator int() { return fIndex; }
    int operator()() { return fIndex; }

    LoopIndex operator=(const LoopIndex &t)
    {
        fLoopSize = t.fLoopSize;
        fIndex = t.fIndex % fLoopSize;
    }

    int operator=(int index)
    {
        if (fLoopSize == 0)
        {
            fIndex = 0;
            return index;
        }
        fIndex = index % fLoopSize;
        return index;
    }

    int operator++(int)
    {
        if (fIndex + 1 < fLoopSize)
        {
            return fIndex++;
        }
        else
        {
            fIndex = 0;
            return fLoopSize - 1;
        }
    }

    int operator--(int)
    {
        if (fIndex - 1 > 0)
        {
            return fIndex--;
        }
        else
        {
            fIndex = fLoopSize - 1;
            return 0;
        }
    }
};

bool operator==(const LoopIndex &t1, const LoopIndex &t2)
{
    if ((t1.fLoopSize == t2.fLoopSize) && (t1.fIndex % t1.fLoopSize == t2.fIndex % t2.fLoopSize))
        return true;
    return false;
}
bool operator!=(const LoopIndex &t1, const LoopIndex &t2) { return !(t1 == t2); }

// const int

// Generate Chain from Given folder
TChain *Generate_Chain(string sFolderName, TChain *Chain);

struct TimeStamp // Time stamp is aligned to board mac0
{
    double fTimeStamp;
    double fCounter;
} gTS;

struct T0DataStat // Get T0Data Statistic information in order to set TS cut condition
{
    double fOffset;
    double fStdDev;
} gT0DataStat;

// Define global TS judgement
bool JudgeMppcTS(int ts0, int ts1)
{
    if (ts0 == 0 || ts1 == 0)
        return false;
    if (TMath::Abs(ts0 - ts1) < 2000) // 800ns
        return true;
    return false;
}

bool JudgeT0TS(double t0ts, int mppcts)
{
    if (t0ts == 0 || mppcts == 0)
        return false;
    if (TMath::Abs((t0ts - gT0DataStat.fOffset) / 1e3 - mppcts) < 3 * gT0DataStat.fStdDev) // about 2500ns
        return true;
    return false;
}

// mppc mac board map
int Mac2Index(UChar_t mac5)
{
    switch (mac5)
    {
    case 0:
        return 0;
    case 1:
        return 1;
    case 2:
        return 2;
    case 5:
        return 3;
    default:
        return -1;
    }
    return -1;
}

UChar_t Index2Mac(int index)
{
    switch (index)
    {
    case 0:
        return 0;
    case 1:
        return 1;
    case 2:
        return 2;
    case 3:
        return 5;
    default:
        return -1;
    }
    return -1;
}

class MppcData // mppc data manager
{
public:
    MppcData(string sFolderName = "./", TTree *tree = 0); // Construct function
    virtual ~MppcData();
    virtual Int_t Cut(Long64_t entry);
    virtual Int_t GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void Init(TTree *tree);
    virtual Bool_t Notify();
    virtual void Show(Long64_t entry = -1);
    virtual int GetEntries() { return fChain->GetEntries(); }

    // Add different versions of Loop function
    virtual void Loop();

    // private:
    int fReadCounter; //!
    int fEntries;     //!

    // Declaration of leaf types
    UChar_t mac5;
    UShort_t chg[32];
    UInt_t ts0;
    UInt_t ts1;
    UInt_t ts0_ref;
    UInt_t ts1_ref;

    // List of branches
    TBranch *b_mac5;    //!
    TBranch *b_chg;     //!
    TBranch *b_ts0;     //!
    TBranch *b_ts1;     //!
    TBranch *b_ts0_ref; //!
    TBranch *b_ts1_ref; //!

    TTree *fChain;  //!pointer to the analyzed TTree or TChain
    Int_t fCurrent; //!current Tree number in a TChain
};

int main(int argc, char **argv)
{
    string sFolder = "./";
    if (argc > 1)
    {
        sFolder = argv[1];
    }

    // TODO: Add Judge whether hfile exists or h is normal
    // Get Stat infomation
    auto hFile = new TFile(Form("%s/Compare.root", sFolder.c_str()));
    auto h = (TH1 *)hFile->Get("h");
    h->Fit("gaus");
    auto f = h->GetFunction("gaus");
    gT0DataStat.fOffset = f->GetParameter(1);
    gT0DataStat.fStdDev = f->GetParameter(2);
    hFile->Close();

    // Read mppc data
    MppcData mppcMG(sFolder);

    // Read t0data
    auto t0File = new TFile(Form("%s/t0data.root", sFolder.c_str()));
    auto t0Tree = (TTree *)t0File->Get("t0data");
    int t0serial;
    double ch5[4];
    t0Tree->SetBranchAddress("serial", &t0serial);
    t0Tree->SetBranchAddress("ch5", &ch5);

    // Prepare match key file
    auto saveFile = new TFile(Form("%s/MatchKey.root", sFolder.c_str()), "recreate");
    auto saveTree = new TTree("key", "key");
    int mppcCount;
    int mppcKey[4], t0Key;
    int mppcMac[4];
    bool t0Flag;
    double ts;
    saveTree->Branch("ts", &ts, "ts/D");
    saveTree->Branch("mppcCount", &mppcCount, "mppcCount/I");
    saveTree->Branch("mppcKey", &mppcKey, "mppcKey[mppcCount]/I");
    saveTree->Branch("mppcMac", &mppcMac, "mppcMac[mppcCount]/I");
    saveTree->Branch("t0Key", &t0Key, "t0Key/I");
    saveTree->Branch("t0Flag", &t0Flag, "t0Flag/B");

    // Start Match
    int singleMppcDataCount = 0; // count mppc data in the loop
    int singleT0DataCount = 0;   // count T0 data in the loop

    int matchedT0DataCount = 0;   // count matched T0Data in matched loop
    int matchedMppcDataCount = 0; // count matched mppc data in matched loop, (no T0 data)

    LoopIndex singleHead(gLoopSize);    // Change only when insert new data
    LoopIndex singleTail(gLoopSize);    // Change only when clear old data
    int singleMppcIndexLoop[gLoopSize]; // Single mppc Read index loop
    int singleT0IndexLoop[gLoopSize];   // Single T0 Read index loop
    double singleTSLoop[gLoopSize];     // single time stamp loop
    DataType singleTypeLoop[gLoopSize]; // single data type loop
    int singleMacLoop[gLoopSize];       // single data Mac loop
    bool singleFlagLoop[gLoopSize];     // single flag loop

    LoopIndex matchedHead(gLoopSize);   // Change only when insert new data
    LoopIndex matchedTail(gLoopSize);   // Change only when clear old data
    double matchedTSLoop[gLoopSize];    // Loop for matched TS
    double t0IndexLoop[gLoopSize];      // Loop for t0data index
    bool t0FlagLoop[gLoopSize];         // Loop for t0 data flag
    double mppcIndexLoop[gLoopSize][4]; // Loop for mppc data index
    int mppcCountLoop[gLoopSize];       // Loop for mppc data counter
    int mppcMacLoop[gLoopSize][4];      // Loop for mppc data Mac

    // Initiate all Loops
    singleHead = 0;
    singleTail = 0;
    matchedHead = 0;
    matchedTail = 0;
    for (int i = 0; i < gLoopSize; i++)
    {
        singleMppcIndexLoop[i] = -1;
        singleT0IndexLoop[i] = -1;
        singleTSLoop[i] = -1;
        singleTypeLoop[i] = undefined;
        singleMacLoop[i] = -1;
        singleFlagLoop[i] = false;

        matchedTSLoop[i] = -1;
        t0IndexLoop[i] = -1;
        for (int j = 0; j < 4; j++)
        {
            mppcIndexLoop[i][j] = -1;
            mppcMacLoop[i][j] = -1;
        }
        t0FlagLoop[i] = false;
        mppcCountLoop[i] = 0;
    }

    // Set start and end index for both data stream
    int mppcStart = 0, mppcEnd = mppcMG.GetEntries();
    int t0Start = 0, t0End = t0Tree->GetEntries();
    int mppctreeIndex = 0, t0treeIndex = 0;

    int mppcReadIndex = 0, t0ReadIndex = 0;

    int loopCount = 0;
    while (1)
    {
        loopCount++;
        int MatchInMatchedLoop = 0;
        int MatchInSingleLoop = 0;

        std::cout << "Loop: " << loopCount << '\t';

        // Compare mppc data count and t0 data count inside single loop
        // if (1)
        if (singleMppcDataCount + matchedMppcDataCount <= singleT0DataCount + matchedT0DataCount)
        {
            std::cout << "Reading mppc data" << std::endl;
            for (int i = 0; i < gLoopSize / 4; i++, mppcReadIndex++)
            {
                if (mppcReadIndex == mppcEnd)
                {
                    break;
                }
                mppcMG.GetEntry(mppcReadIndex);
                int mppcts = mppcMG.ts0;
                if (mppcts == 0)
                    continue;

                // Search in matched Loop
                bool matchFlag = 0;

                if (matchedTail != matchedHead)
                {
                    int i = 0;
                    for (LoopIndex loopidx = matchedTail; loopidx != matchedHead + 1 || i == 0; loopidx++, i++)
                    {
                        int arrayidx = loopidx();
                        if (!t0FlagLoop[arrayidx] && mppcCountLoop[arrayidx] == 0)
                            continue;

                        if (JudgeMppcTS(mppcts, matchedTSLoop[arrayidx])) // in the case of matched
                        {
                            int boardNo = Mac2Index(mppcMG.mac5);
                            int boardCount = mppcCountLoop[arrayidx];
                            if (boardCount < 4)
                            {
                                // Check all mac5 in this matched data
                                bool mac5CheckFlag = 1;
                                for (int i = 0; i < boardCount; i++)
                                {
                                    UChar_t mac5 = mppcMacLoop[arrayidx][i];
                                    if (mppcMG.mac5 == mac5)
                                    {
                                        mac5CheckFlag = 0;
                                    }
                                }
                                if (!mac5CheckFlag)
                                    continue;

                                matchedTSLoop[arrayidx] = (matchedTSLoop[arrayidx] * boardCount + mppcts) / (boardCount + 1);

                                mppcMacLoop[arrayidx][boardCount] = mppcMG.mac5;
                                mppcIndexLoop[arrayidx][mppcCountLoop[arrayidx]++] = mppcReadIndex;

                                matchFlag = 1;
                                break;
                            }
                        }
                    }
                }

                if (matchFlag)
                    continue;

                // If not found in matched loop, serach in single loop
                if (singleHead != singleTail)
                {
                    int i = 0;
                    for (LoopIndex loopidx = singleTail; loopidx != singleHead + 1 || i == 0; loopidx++, i++)
                    {
                        int arrayidx = loopidx();
                        if (!singleFlagLoop[arrayidx])
                            continue;

                        auto datatype = singleTypeLoop[arrayidx];
                        if (datatype == undefined)
                            continue;

                        bool tsFlag;
                        if (datatype == mppc)
                        {
                            tsFlag = JudgeMppcTS(singleTSLoop[arrayidx], mppcMG.ts0);
                        }
                        else if (datatype == t0)
                        {
                            tsFlag = JudgeT0TS(singleTSLoop[arrayidx], mppcMG.ts0);
                        }

                        if (!tsFlag)
                            continue;

                        // in the case of matched
                        matchFlag = true;
                        // before insert into matched loop, save and clear replaced element
                        if ((matchedHead + 1) == matchedTail) // if head+1==tail, shows that the loop is full
                        {
                            // Save matched data
                            int saveidx = matchedTail();
                            ts = matchedTSLoop[saveidx];
                            t0Key = t0IndexLoop[saveidx];
                            t0Flag = t0FlagLoop[saveidx];

                            mppcCount = mppcCountLoop[saveidx];
                            for (int j = 0; j < mppcCount; j++)
                            {
                                mppcKey[j] = mppcIndexLoop[saveidx][j];
                                mppcMac[j] = mppcMacLoop[saveidx][j];
                            }
                            saveTree->Fill();

                            // Clear matched data
                            if (t0FlagLoop[saveidx])
                            {
                                matchedT0DataCount--;
                            }
                            else
                            {
                                matchedMppcDataCount--;
                            }

                            matchedTSLoop[saveidx] = -1;
                            t0IndexLoop[saveidx] = -1;
                            t0FlagLoop[saveidx] = false;
                            mppcCountLoop[saveidx] = 0;
                            for (int j = 0; j < 4; j++)
                            {
                                mppcIndexLoop[saveidx][j] = -1;
                                mppcMacLoop[saveidx][j] = -1;
                            }

                            matchedTail++;
                        }

                        int insertidx;
                        if (matchedHead == matchedTail) // if head == tail, shows that this loop is empty
                        {
                            insertidx = matchedHead();
                            matchedHead++;
                        }
                        else
                        {
                            matchedHead++;
                            insertidx = matchedHead();
                        }

                        if (datatype == mppc)
                        {
                            double matchedTS = (singleTSLoop[arrayidx] + mppcMG.ts0) / 2.0;
                            matchedTSLoop[insertidx] = matchedTS;
                            mppcCountLoop[insertidx] = 2;
                            mppcIndexLoop[insertidx][0] = singleMppcIndexLoop[arrayidx];
                            mppcMacLoop[insertidx][0] = singleMacLoop[arrayidx];
                            mppcIndexLoop[insertidx][1] = mppcReadIndex;
                            mppcMacLoop[insertidx][1] = mppcMG.mac5;

                            t0FlagLoop[insertidx] = false;
                            t0IndexLoop[insertidx] = -1;

                            singleMppcDataCount--;
                            matchedMppcDataCount++;
                        }
                        else if (datatype == t0)
                        {
                            matchedTSLoop[insertidx] = mppcMG.ts0;
                            mppcCountLoop[insertidx] = 1;
                            mppcIndexLoop[insertidx][0] = mppcReadIndex;
                            mppcMacLoop[insertidx][0] = mppcMG.mac5;

                            t0FlagLoop[insertidx] = true;
                            t0IndexLoop[insertidx] = singleT0IndexLoop[arrayidx];

                            singleT0DataCount--;
                            matchedT0DataCount++;
                        }

                        // Clear this element in single loop
                        singleMppcIndexLoop[arrayidx] = -1;
                        singleT0IndexLoop[arrayidx] = -1;
                        singleTSLoop[arrayidx] = -1;
                        singleTypeLoop[arrayidx] = undefined;
                        singleMacLoop[arrayidx] = -1;
                        singleFlagLoop[arrayidx] = false;

                        break;
                    }
                }

                if (matchFlag)
                    continue;

                // If still not found, insert into single loop
                if ((singleHead + 1) == singleTail) // in the case of loop is full
                {
                    // throw warning:
                    int saveidx = singleTail();

                    if (singleFlagLoop[saveidx] == true)
                    // if (0)
                    {

                        std::cerr << "Warning: This data is not matched with any data" << std::endl;
                        std::cerr << "Datatype: " << singleTypeLoop[saveidx] << std::endl;
                        std::cerr << "Index: " << ((singleTypeLoop[saveidx] == mppc) ? singleMppcIndexLoop[saveidx] : singleT0IndexLoop[saveidx]) << std::endl;
                        std::cerr << "mac5: " << singleMacLoop[saveidx] << std::endl;
                        std::cerr << "TS: " << singleTSLoop[saveidx] << std::endl;
                    }

                    // Clear
                    if (singleTypeLoop[saveidx] == mppc && singleFlagLoop[saveidx] == 1)
                    {
                        singleMppcDataCount--;
                    }
                    else if (singleTypeLoop[saveidx] == t0 && singleFlagLoop[saveidx] == 1)
                    {
                        singleT0DataCount--;
                    }
                    singleMppcIndexLoop[saveidx] = -1;
                    singleT0IndexLoop[saveidx] = -1;
                    singleTSLoop[saveidx] = -1;
                    singleTypeLoop[saveidx] = undefined;
                    singleMacLoop[saveidx] = -1;
                    singleFlagLoop[saveidx] = false;

                    singleTail++;
                }

                int insertidx;
                if (singleHead == singleTail)
                {
                    insertidx = singleHead();
                    singleHead++;
                }
                else
                {
                    singleHead++;
                    insertidx = singleHead();
                }

                singleMppcIndexLoop[insertidx] = mppcReadIndex;
                singleTSLoop[insertidx] = mppcMG.ts0;
                singleTypeLoop[insertidx] = mppc;
                singleMacLoop[insertidx] = mppcMG.mac5;
                singleFlagLoop[insertidx] = true;

                singleMppcDataCount++;
            }
        }
        else
        // if (1)
        {
            std::cout << "Reading t0 data" << std::endl;
            for (int i = 0; i < gLoopSize / 4; i++, t0ReadIndex++)
            {
                if (t0ReadIndex == t0End)
                {
                    break;
                }
                t0Tree->GetEntry(t0ReadIndex);

                double t0TS = ch5[0];

                // Search in matched Loop
                bool matchFlag = 0;

                if (matchedTail != matchedHead)
                {
                    int i = 0;
                    for (LoopIndex loopidx = matchedTail; loopidx != matchedHead + 1 || i == 0; loopidx++, i++)
                    {
                        int arrayidx = loopidx();
                        if (!t0FlagLoop[arrayidx] && mppcCountLoop[arrayidx] == 0)
                            continue;
                        // cout << t0TS << '\t' << matchedTSLoop[arrayidx] << endl;
                        if (JudgeT0TS(t0TS, matchedTSLoop[arrayidx])) // in the case of matched
                        {
                            if (t0FlagLoop[arrayidx])
                                continue;

                            int boardCount = mppcCountLoop[arrayidx];
                            t0IndexLoop[arrayidx] = t0ReadIndex;
                            t0FlagLoop[arrayidx] = true;

                            matchFlag = 1;

                            matchedMppcDataCount--;
                            matchedT0DataCount++;

                            // cout << "Matched: " << endl;
                            break;
                        }
                    }
                }
                if (matchFlag)
                    continue;

                // If not found in matched loop, serach in single loop
                if (singleHead != singleTail)
                {
                    int i = 0;
                    for (LoopIndex loopidx = singleTail; loopidx != singleHead + 1 || i == 0; loopidx++, i++)
                    {
                        int arrayidx = loopidx();
                        if (!singleFlagLoop[arrayidx])
                            continue;

                        auto datatype = singleTypeLoop[arrayidx];
                        if (datatype == undefined || datatype == t0)
                            continue;

                        bool tsFlag;
                        if (datatype == mppc)
                        {
                            tsFlag = JudgeT0TS(t0TS, singleTSLoop[arrayidx]);
                        }

                        if (!tsFlag)
                            continue;

                        // in the case of matched
                        // before insert into matched loop, save and clear replaced element
                        if (matchedTail == (matchedHead + 1)) // if head+1==tail, shows that the loop is full
                        {
                            // Save matched data
                            int saveidx = matchedTail();
                            ts = matchedTSLoop[saveidx];
                            t0Key = t0IndexLoop[saveidx];
                            t0Flag = t0FlagLoop[saveidx];

                            mppcCount = mppcCountLoop[saveidx];
                            for (int j = 0; j < mppcCount; j++)
                            {
                                mppcKey[j] = mppcIndexLoop[saveidx][j];
                                mppcMac[j] = mppcMacLoop[saveidx][j];
                            }
                            saveTree->Fill();

                            // Clear matched data
                            if (t0FlagLoop[saveidx] == false)
                            {
                                matchedMppcDataCount--;
                            }
                            else
                            {
                                matchedT0DataCount--;
                            }

                            matchedTSLoop[saveidx] = -1;
                            t0IndexLoop[saveidx] = -1;
                            t0FlagLoop[saveidx] = false;
                            mppcCountLoop[saveidx] = 0;
                            for (int j = 0; j < 4; j++)
                            {
                                mppcIndexLoop[saveidx][j] = -1;
                                mppcMacLoop[saveidx][j] = -1;
                            }

                            matchedTail++;
                        }

                        int insertidx;
                        if (matchedHead == matchedTail) // if head == tail, shows that this loop is empty
                        {
                            insertidx = matchedHead();
                            matchedHead++;
                        }
                        else
                        {
                            matchedHead++;
                            insertidx = matchedHead();
                        }

                        matchedTSLoop[insertidx] = singleTSLoop[arrayidx];
                        mppcCountLoop[insertidx] = 1;
                        mppcIndexLoop[insertidx][0] = singleMppcIndexLoop[arrayidx];
                        mppcMacLoop[insertidx][0] = singleMacLoop[arrayidx];

                        t0FlagLoop[insertidx] = true;
                        t0IndexLoop[insertidx] = t0ReadIndex;

                        singleMppcDataCount--;
                        matchedT0DataCount++;

                        // Clear this element in single loop
                        singleMppcIndexLoop[arrayidx] = -1;
                        singleT0IndexLoop[arrayidx] = -1;
                        singleTSLoop[arrayidx] = -1;
                        singleTypeLoop[arrayidx] = undefined;
                        singleMacLoop[arrayidx] = -1;
                        singleFlagLoop[arrayidx] = false;

                        matchFlag = true;
                        break;
                    }
                }
                if (matchFlag)
                    continue;

                // If still not found, insert into single loop
                if (singleTail == (singleHead + 1)) // in the case of loop is full
                {
                    // throw warning:
                    int saveidx = singleTail();
                    if (singleTypeLoop[saveidx] != undefined)
                    {
                        std::cout << "Warning: This data is not matched with any data" << std::endl;
                        std::cout << "Datatype: " << singleTypeLoop[saveidx] << std::endl;
                        std::cout << "Index: " << ((singleTypeLoop[saveidx] == mppc) ? singleMppcIndexLoop[saveidx] : singleT0IndexLoop[saveidx]) << std::endl;
                        std::cout << "mac5: " << singleMacLoop[saveidx] << std::endl;
                        std::cout << "TS: " << singleTSLoop[saveidx] << std::endl;
                    }

                    // Clear
                    singleMppcIndexLoop[saveidx] = -1;
                    singleT0IndexLoop[saveidx] = -1;
                    singleTSLoop[saveidx] = -1;
                    singleTypeLoop[saveidx] = undefined;
                    singleMacLoop[saveidx] = -1;
                    singleFlagLoop[saveidx] = false;

                    if (singleTypeLoop[saveidx] == mppc && singleFlagLoop[saveidx] == 1)
                    {
                        singleMppcDataCount--;
                    }
                    if (singleTypeLoop[saveidx] == t0 && singleFlagLoop[saveidx] == 1)
                    {
                        singleT0DataCount--;
                    }

                    singleTail++;
                }

                int insertidx;
                if (singleHead == singleTail)
                {
                    insertidx = singleHead();
                    singleHead++;
                }
                else
                {
                    singleHead++;
                    insertidx = singleHead();
                }

                singleT0IndexLoop[insertidx] = t0ReadIndex;
                singleTSLoop[insertidx] = t0TS;
                singleTypeLoop[insertidx] = t0;
                singleFlagLoop[insertidx] = true;

                singleT0DataCount++;
            }
        }

        {
            std::cout << "mppc: " << mppcReadIndex << '\t' << "t0: " << t0ReadIndex << std::endl;
            std::cout << "Mppc in single loop: " << singleMppcDataCount << "\t T0 in single loop: " << singleT0DataCount << std::endl;
            std::cout << "Mppc in matched loop: " << matchedMppcDataCount << "\t T0 in matched loop: " << matchedT0DataCount << std::endl
                      << std::endl;
        }
        // if (loopCount == 10000)
        //     return 0;

        if (mppcReadIndex == mppcEnd || t0ReadIndex == t0End)
        {
            for (LoopIndex loopidx = matchedTail; loopidx != matchedHead; loopidx++)
            {
                // Save matched data
                int saveidx = loopidx();
                ts = matchedTSLoop[saveidx];
                t0Key = t0IndexLoop[saveidx];
                t0Flag = t0FlagLoop[saveidx];

                mppcCount = mppcCountLoop[saveidx];
                for (int j = 0; j < mppcCount; j++)
                {
                    mppcKey[j] = mppcIndexLoop[saveidx][j];
                    mppcMac[j] = mppcMacLoop[saveidx][j];
                }
                saveTree->Fill();
            }
            break;
        }
    }

    saveTree->Write();
    saveFile->Close();
}

TChain *Generate_Chain(string sFolderName, TChain *Chain)
{
    gSystem->Exec(Form("ls %s/mppc*.root > .~filelist", sFolderName.c_str()));

    auto ch = Chain;
    if (ch == NULL)
    {
        ch = new TChain("mppc");
    }

    ifstream file_list(".~filelist");
    for (int i = 0; file_list.is_open() && file_list.eof() == false; i++)
    {
        string s_temp;
        file_list >> s_temp;
        if (s_temp.find(".root") == string::npos)
        {
            continue;
        }
        else if (s_temp.find("mppc") == string::npos)
        {
            continue;
        }
        else if (s_temp.find("Histo") != string::npos)
        {
            continue;
        }

        cout << "File: " << s_temp << " Read" << endl;
        ch->Add(s_temp.c_str());
    }
    cout << "Totally get " << ch->GetEntries() << " Entries" << endl;
    gSystem->Exec("rm .~filelist");
    return ch;
}

MppcData::MppcData(string sFolderName, TTree *tree) : fChain(0)
{

    // Construct fBoardNum & Reset clock counter;
    fReadCounter = 0;
    fEntries = 0;
    if (tree == 0)
    {
        auto chain = new TChain("mppc");
        Generate_Chain(sFolderName, chain);
        tree = chain;
    }
    if (tree)
        Init(tree);
    if (b_mac5)
        fEntries = fChain->GetEntries();
}

MppcData::~MppcData()
{
    if (!fChain)
        return;
    delete fChain->GetCurrentFile();
}

Int_t MppcData::GetEntry(Long64_t entry)
{
    if (!fChain)
        return 0;
    return fChain->GetEntry(entry);
}
Long64_t MppcData::LoadTree(Long64_t entry)
{
    if (!fChain)
        return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0)
        return centry;
    if (fChain->GetTreeNumber() != fCurrent)
    {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void MppcData::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree)
        return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("mac5", &mac5, &b_mac5);
    fChain->SetBranchAddress("chg", chg, &b_chg);
    fChain->SetBranchAddress("ts0", &ts0, &b_ts0);
    fChain->SetBranchAddress("ts1", &ts1, &b_ts1);
    fChain->SetBranchAddress("ts0_ref", &ts0_ref, &b_ts0_ref);
    fChain->SetBranchAddress("ts1_ref", &ts1_ref, &b_ts1_ref);
    Notify();
}

Bool_t MppcData::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void MppcData::Show(Long64_t entry)
{
    if (!fChain)
        return;
    fChain->Show(entry);
}
Int_t MppcData::Cut(Long64_t entry)
{
    if ((ts0 == 0) || (ts1 == 0))
    {
        return -1;
    }
    return 1;
}

void MppcData::Loop()
{
    if (fChain == 0)
        return;
}
