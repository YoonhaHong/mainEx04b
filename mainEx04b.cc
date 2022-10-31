#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "TFile.h"
#include "Pythia8Plugins/FastJet3.h"
#define RMAX 6
#define MOTHER 30
using namespace Pythia8;

int main(int argc, char *argv[]) {

	if ( argc<2 ){
		cout << "Usage: ./mainEx04b cut_jetpt" << endl;
		return -1;
	}

	float cut_jetpt = atof(argv[1]);
	cout << "Cut Jet pT: " << cut_jetpt << endl;

	Pythia pythia;

	Event& event = pythia.event;

	pythia.readFile("mainEx04b.cfg");

	int nEvent = pythia.mode("Main:numberOfEvents");
	int nAbort = pythia.mode("Main:timesAllowErrors");

	pythia.init();
	fastjet::JetDefinition *jetDef[RMAX];
	for(int R=0; R<RMAX; R++){jetDef[R]= new fastjet::JetDefinition(fastjet::genkt_algorithm, (double)(R+2)/(double)10.0, -1);} //0.2~0.7

    std::vector <fastjet::PseudoJet> fjInputs;
	std::vector <fastjet::PseudoJet> fjInputsChg;
	const float const_pt_cut = 0.15;
    const float const_eta_cut = 4;

	int i_np;
	int i_p_id[5000]; 
    int i_p_status[5000];
	bool b_p_chg[5000];
	float f_p_pt[5000];
	float f_p_eta[5000];
	float f_p_phi[5000];
	float f_p_e[5000];
	float f_p_m[5000];
	float f_p_t[5000];

	int i_nchgjet[RMAX];
	float f_chgjet_pt[RMAX][20];
	float f_chgjet_eta[RMAX][20];
	float f_chgjet_phi[RMAX][20];

    int i_nijet[RMAX];
    float f_ijet_pt[RMAX][20];
    float f_ijet_eta[RMAX][20];
    float f_ijet_phi[RMAX][20];
    
	auto T = new TTree("T","Pythia event");

	T->Branch("np",&i_np,"np/I");
	T->Branch("p_id",i_p_id,"p_id[np]/I");
    T->Branch("p_status",i_p_status,"p_status[np]/I");
	T->Branch("p_chg",b_p_chg,"p_chg[np]/O");
	T->Branch("p_pt",f_p_pt,"p_pt[np]/F");
	T->Branch("p_eta",f_p_eta,"p_eta[np]/F");
	T->Branch("p_phi",f_p_phi,"p_phi[np]/F");
	T->Branch("p_e",f_p_e,"p_e[np]/F");
	T->Branch("p_m",f_p_m,"p_m[np]/F");
	T->Branch("p_t",f_p_t,"p_t[np]/F");

    
	TString chgjet;
	for(int R=0; R<RMAX; R++){
		chgjet.Form("chgjetR%02d",R+2);
		T->Branch("n"+chgjet,&i_nchgjet[R],"n"+chgjet+"/I");
		T->Branch(chgjet+"_pt",f_chgjet_pt[R],chgjet+"_pt[n"+chgjet+"]/F");
		T->Branch(chgjet+"_eta",f_chgjet_eta[R],chgjet+"_eta[n"+chgjet+"]/F");
		T->Branch(chgjet+"_phi",f_chgjet_phi[R],chgjet+"_phi[n"+chgjet+"]/F");
	}
	TString ijet;
	for(int R=0; R<RMAX; R++){
		ijet.Form("ijetR%02d",R+2);
		T->Branch("n"+ijet,&i_nijet[R],"n"+ijet+"/I");
		T->Branch(ijet+"_pt",f_ijet_pt[R],ijet+"_pt[n"+ijet+"]/F");
		T->Branch(ijet+"_eta",f_ijet_eta[R],ijet+"_eta[n"+ijet+"]/F");
		T->Branch(ijet+"_phi",f_ijet_phi[R],ijet+"_phi[n"+ijet+"]/F");
	}

  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

		if (!pythia.next()) {
			if (++iAbort < nAbort) continue;
			cout << " Event generation aborted prematurely, owing to error!\n";
			break;
		}
		fjInputsChg.resize(0);
        fjInputs.resize(0);

		for (int i = 0; i < event.size(); ++i) {
			if ( !(event[i].isFinal()) ) continue; 
			if ( fabs(event[i].eta())>const_eta_cut ) continue;

			fastjet::PseudoJet particleTemp = event[i];

            fjInputs.push_back(particleTemp);
            if ( !(event[i].isCharged()) ) continue;
			if ( event[i].pT()<const_pt_cut ) continue;
			fjInputsChg.push_back(particleTemp);
		}

    vector <fastjet::PseudoJet> inclusiveChgJets[RMAX], sortedChgJets[RMAX];
    fastjet::ClusterSequence* clustSeqChg[RMAX];
    vector <fastjet::PseudoJet> inclusiveJets[RMAX], sortedJets[RMAX];
    fastjet::ClusterSequence* clustSeq[RMAX];


		for(int R=0; R<RMAX; R++){
			clustSeqChg[R] = new fastjet::ClusterSequence(fjInputsChg,*jetDef[R]);
			inclusiveChgJets[R] = clustSeqChg[R]->inclusive_jets(cut_jetpt);
			sortedChgJets[R] = sorted_by_pt(inclusiveChgJets[R]);
		
            clustSeq[R] = new fastjet::ClusterSequence(fjInputs,*jetDef[R]);
			inclusiveJets[R] = clustSeq[R]->inclusive_jets(cut_jetpt);
			sortedJets[R] = sorted_by_pt(inclusiveJets[R]);
	}

		for(int R=0; R<RMAX; R++){
			i_nchgjet[R]=0;
            i_nijet[R]=0;
			for (int ii=0; ii<10; ii++){
                f_chgjet_pt[R][ii] = f_chgjet_eta[R][ii] = f_chgjet_phi[R][ii] = -999;
                f_ijet_pt[R][ii] = f_ijet_eta[R][ii] = f_ijet_phi[R][ii] = -999;           }
		}

		int i_nijet_all = 0;
	
		for(int R=0; R<RMAX; R++){
			for (int i = 0; i < int(sortedChgJets[R].size()); ++i) {

				if ( fabs(sortedChgJets[R][i].rap())>3) continue;

				f_chgjet_pt[R][i_nchgjet[R]] = sortedChgJets[R][i].perp();
				f_chgjet_eta[R][i_nchgjet[R]] = sortedChgJets[R][i].rap();
				f_chgjet_phi[R][i_nchgjet[R]] = sortedChgJets[R][i].phi_std();

				i_nchgjet[R]++;
                i_nijet_all++;
			}
			for (int j = 0; j < int(sortedJets[R].size()); ++j) {

				if ( fabs(sortedJets[R][j].rap())>3) continue;

				f_ijet_pt[R][i_nijet[R]] = sortedJets[R][j].perp();
				f_ijet_eta[R][i_nijet[R]] = sortedJets[R][j].rap();
				f_ijet_phi[R][i_nijet[R]] = sortedJets[R][j].phi_std();

				i_nijet[R]++;
                i_nijet_all++;
			}
    	}

		if ( i_nijet_all<1 ) continue;


		i_np = 0;
		for (int i = 0; i < event.size(); ++i) {

            if ( event[i].isFinal() ) continue; //only for tracking mother

			float tmp_pt = event[i].pT();
			float tmp_eta = event[i].eta();
			float tmp_phi = event[i].phi();
            int id = event[i].id();

			if ( fabs(tmp_eta)>const_eta_cut ) continue;
			if ( abs(id)==12 ||abs(id)==14 || abs(id)==16 ) continue;
			i_p_id[i_np] =  id;
            i_p_status[i_np] = event[i].status();
			b_p_chg[i_np] = event[i].isCharged();
			f_p_pt[i_np] = tmp_pt;
			f_p_eta[i_np] = tmp_eta;
			f_p_phi[i_np] = tmp_phi;
			f_p_e[i_np] = event[i].e();
			f_p_m[i_np] = event[i].m();
			f_p_t[i_np] = event[i].tProd();
            
			i_np++;
		}//i

		T->Fill();
		for(int R=0; R<RMAX; R++){
			delete clustSeqChg[R];
			delete clustSeq[R];
		}  
  }

  // Final statistics. Normalize and output histograms.
  pythia.stat();

	TFile *outfile = new TFile("Pythia8_event.root","recreate");
	T->Write();
	outfile->Close();

  // Done.
  return 0;

}
