#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"
#include "TCanvas.h"

using namespace std;

struct T0Data
{
    int serial;
    double data[4];
};

ifstream &ReadT0Line(ifstream &fin, vector<double> &t0data)
{
    double tl, th, pl, ph;
    char c;

    int serial;
    fin >> serial;
    fin >> c;

    fin >> tl;
    fin >> c;
    fin >> th;
    fin >> c;
    fin >> pl;
    fin >> c;
    fin >> ph;

    t0data.push_back(tl);
    return fin;
}

ifstream &ReadT0Line(ifstream &fin, T0Data &t0data)
{
    double tl, th, pl, ph;
    char c;

    int serial;
    fin >> serial;
    fin >> c;

    fin >> tl;
    fin >> c;
    fin >> th;
    fin >> c;
    fin >> pl;
    fin >> c;
    fin >> ph;

    t0data.serial = serial;
    t0data.data[0] = tl;
    t0data.data[1] = th;
    t0data.data[2] = pl;
    t0data.data[3] = ph;
    return fin;
}

bool JudgeReadable(ifstream &fin)
{
    return fin.good() && fin.is_open() && !fin.eof();
}

// ifstream &ReadT0Line(ifstream &fin, vector<double> &t0data)
// {
//     double tl, th, pl, ph;
//     char c;

//     // int serial;
//     // fin >> serial;
//     // fin >> c;

//     fin >> tl;

//     t0data.push_back(tl);
//     return fin;
// }

void Compare(string sFolderName = "./", string sDate = "20210119")
{
    auto h = new TH1D("h", "h", 200, -2.5e3, 0);

    auto file = new TFile((sFolderName + "mppc.root").c_str());
    auto tree = (TTree *)file->Get("mppc");

    UChar_t mac5;
    UInt_t ts;
    tree->SetBranchAddress("mac5", &mac5);
    tree->SetBranchAddress("ts0", &ts);

    vector<double> mppcdata;
    vector<int> mppcSerial;
    for (int i = 0; i < tree->GetEntries(); i++)
    {
        tree->GetEntry(i);
        if (mac5 == 0 && ts != 0)
        {
            mppcdata.push_back(ts);
            mppcSerial.push_back(i);
        }
    }
    file->Close();

    vector<double> t0data;
    ifstream fin(Form("%s/CH5_Data_%s_correct_calculated.txt", sFolderName.c_str(), sDate.c_str()));
    for (; fin.good() && !fin.eof() && fin.is_open();)
    {
        ReadT0Line(fin, t0data);
    }
    fin.close();

    int Count = 0;
    cout << "mppc size: " << mppcdata.size() << '\t' << "t0 size: " << t0data.size() << endl;

    auto saveFile = new TFile(Form("%s/Compare.root", sFolderName.c_str()), "recreate");
    auto saveTree = new TTree("key", "key");
    int mppcKey;
    int t0Key;
    saveTree->Branch("mppcKey", &mppcKey, "mppcKey/I");
    saveTree->Branch("t0Key", &t0Key, "t0Key/I");

    for (int i = 0; i < mppcdata.size(); i++)
    {
        if (i % 100 == 0)
        {
            cout << "Processing: " << i << endl;
        }
        double tsmppc = mppcdata[i];
        for (int j = 0; j < t0data.size(); j++)
        {
            double tst0 = t0data[j];

            double temp = tsmppc - tst0 / 1e3;
            if ((temp <= 0 && temp > -2.5e3))
            {
                if (tst0 == 0)
                    continue;

                mppcKey = mppcSerial[i];
                t0Key = j;
                saveTree->Fill();
                h->Fill(tsmppc - tst0 / 1e3);
                cout << Count++ << '\t' << i << '\t' << tsmppc << '\t' << j << '\t' << tst0 << '\t' << temp << endl;
            }
        }
    }

    saveFile->cd();
    saveTree->Write();
    h->Write();
    auto c = new TCanvas("c", "c", 1);
    h->Draw();
    c->SaveAs("Compare.jpg");

    saveFile->Close();
}

int main(int argc, char **argv)
{
    string sFolder = "./";
    if (argc > 1)
    {
        sFolder = argv[1];
    }

    string date = "20210119";
    if (argc > 2)
    {
        date = argv[2];
    }

    Compare(sFolder, date);
    // Init ifstream
    ifstream fin1(Form("%s/CH%d_Data_%s_correct_calculated.txt", sFolder.c_str(), 1, date.c_str()));
    ifstream fin2(Form("%s/CH%d_Data_%s_correct_calculated.txt", sFolder.c_str(), 2, date.c_str()));
    ifstream fin3(Form("%s/CH%d_Data_%s_correct_calculated.txt", sFolder.c_str(), 3, date.c_str()));
    ifstream fin4(Form("%s/CH%d_Data_%s_correct_calculated.txt", sFolder.c_str(), 4, date.c_str()));
    ifstream fin5(Form("%s/CH%d_Data_%s_correct_calculated.txt", sFolder.c_str(), 5, date.c_str()));

    // Trans to root file
    auto file = new TFile(Form("%s/t0data.root", sFolder.c_str()), "recreate");
    auto tree = new TTree("t0data", "t0data");
    int serial;
    int serial1, serial2, serial3, serial4, serial5;
    double ch1[4], ch2[4], ch3[4], ch4[4], ch5[4];
    tree->Branch("serial", &serial, "serial/I");
    tree->Branch("ch1", &ch1, "ch1[4]/D");
    tree->Branch("ch2", &ch2, "ch2[4]/D");
    tree->Branch("ch3", &ch3, "ch3[4]/D");
    tree->Branch("ch4", &ch4, "ch4[4]/D");
    tree->Branch("ch5", &ch5, "ch5[4]/D");

    T0Data t0temp;
    for (;;)
    {
        ReadT0Line(fin1, t0temp);
        serial1 = t0temp.serial;
        for (int i = 0; i < 4; i++)
        {
            ch1[i] = t0temp.data[i];
        }
        if (!JudgeReadable(fin1))
            break;

        ReadT0Line(fin2, t0temp);
        serial2 = t0temp.serial;
        for (int i = 0; i < 4; i++)
        {
            ch2[i] = t0temp.data[i];
        }
        if (!JudgeReadable(fin2))
            break;

        ReadT0Line(fin3, t0temp);
        serial3 = t0temp.serial;
        for (int i = 0; i < 4; i++)
        {
            ch3[i] = t0temp.data[i];
        }
        if (!JudgeReadable(fin3))
            break;

        ReadT0Line(fin4, t0temp);
        serial4 = t0temp.serial;
        for (int i = 0; i < 4; i++)
        {
            ch4[i] = t0temp.data[i];
        }
        if (!JudgeReadable(fin4))
            break;

        ReadT0Line(fin5, t0temp);
        serial5 = t0temp.serial;
        for (int i = 0; i < 4; i++)
        {
            ch5[i] = t0temp.data[i];
        }
        if (!JudgeReadable(fin5))
            break;

        if (serial1 == serial2 && serial2 == serial3 && serial3 == serial4 && serial4 == serial5)
        {
            serial = serial1;
            tree->Fill();
        }
        else
        {
            cout << "Warning: serial error" << serial1 << '\t' << serial2 << '\t' << serial3 << '\t' << serial4 << '\t' << serial5 << endl;
            break;
        }
    }
    tree->Write();
    file->Close();
}