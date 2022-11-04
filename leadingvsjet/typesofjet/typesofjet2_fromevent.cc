//extended to 0<jet pT<140
//leading particle finding range is fixed : dR<0.7

#include <iostream>
#include <vector>
#include <string>
#define pi 3.141592
using namespace std;
class Pseudo:public TVector3
{
	public:
		Pseudo(float pt=0, float eta=0, float phi=0)
		{
			this->SetPtEtaPhi(pt,eta,phi);
		}
};

inline float getjt(Pseudo jet, Pseudo part) {return (jet.Cross(part)).Mag()/jet.Mag();}

inline float getz(Pseudo jet, Pseudo part) {return jet*part/jet.Mag2();}

float getr(float eta1, float eta2, float phi1, float phi2){ 
	float deta=eta1-eta2;
	float dphi=phi1-phi2;
	if(dphi<-pi) dphi+=2*pi;
	else if(dphi>pi) dphi-=2*pi;
	return sqrt(deta*deta+dphi*dphi);
}
void typesofjet2_fromevent(int Rindex=5, TString inputpath="../HardQCD:all/allfull_MIF.root") 
//5 -> R=0.7 "../../gg2qqbarfull_MIF.root"
{
    //int mode=inputpath.Contains("qq"); //if quark jet mode =1 
    int eventmode = 99;
    if( inputpath.Contains("gg2qqbar") ) eventmode = 1;
    else if( inputpath.Contains("gg2gg") ) eventmode = 2;
    else eventmode = 0;

    float *p_pt=new float[5000];
	float *p_eta=new float[5000];
	float *p_phi=new float[5000];
	bool *p_chg=new bool[5000];
	int *p_id=new int[5000];
    int *p_status=new int[5000];

	float *ijet_pt=new float[20];
	float *ijet_eta=new float[20];
	float *ijet_phi=new float[20];

	int partnum;
	int ijetnum=0;
	int nijet=0,nqtaggedijet=0, ngtaggedijet=0, netaggedijet=0;
    int areaijet[3]={0,0,0}; //numbers of jet pT: 10~40, leading particle pT:20~40

    TH2F *qpt = new TH2F("quarkjet_particle_pt","quarkjet_particle_pt",1400,0,140,1400,0,140);
    TH2F *qeta = new TH2F("quarkjet_particle_eta","quarkjet_particle_eta",60,-3,3,80,-4,4);
    TH2F *qphi = new TH2F("quarkjet_particle_phi","quarkjet_particle_phi",60,-pi,pi,60,-pi,pi);
    TH1F *qdr = new TH1F("quarkdr","quarkdr",70,0.,0.7);

	TH2F *gpt = new TH2F("gluonjet_particle_pt","gluonjet_particle_pt",1400,0,140,1400,0,140);
    TH2F *geta = new TH2F("gluonjet_particle_eta","gluonjet_particle_eta",60,-3,3,80,-4,4);
    TH2F *gphi = new TH2F("gluonjet_particle_phi","gluonjet_particle_phi",60,-pi,pi,60,-pi,pi);
    TH1F *gdr = new TH1F("gluondr","gluondr",70,0.,0.7);

    TH2F *ept = new TH2F("elsejet_particle_pt","elsejet_particle_pt",1400,0,140,1400,0,140);
    TH2F *eeta = new TH2F("elsejet_particle_eta","elsejet_particle_eta",60,-3,3,80,-4,4);
    TH2F *ephi = new TH2F("elsejet_particle_phi","elsejet_particle_phi",60,-pi,pi,60,-pi,pi);
    TH1F *edr = new TH1F("elsedr","elsedr",70,0.,0.7);

    TFile *f = new TFile(inputpath, "read");
	TTree *T = (TTree*)f ->Get("T");
	T->SetBranchAddress("p_pt",p_pt);
	T->SetBranchAddress("p_id",p_id);
    T->SetBranchAddress("p_status",p_status);
	T->SetBranchAddress("np",&partnum);
	T->SetBranchAddress("p_eta",p_eta);
	T->SetBranchAddress("p_phi",p_phi);
	T->SetBranchAddress("p_chg",p_chg);

    TString sijet;
	sijet.Form("ijetR%02d",Rindex+2);
	T->SetBranchAddress(sijet+"_pt",ijet_pt);
    T->SetBranchAddress(sijet+"_eta",ijet_eta);
    T->SetBranchAddress(sijet+"_phi",ijet_phi);
    T->SetBranchAddress("n"+sijet,&ijetnum);

	
    
	float jet_R_size = 0.1*(2.0+(float)Rindex); //R=0.7(Rindex=5) eta=-0.2~0.2
	float leading_pt; 
    int leading_index;
    int tagged_id;
    float dr;

	const int entrysum=T->GetEntries();

	for (int ien=0; ien<entrysum; ien++){
		T->GetEntry(ien);

		if ( ijetnum<1 ) continue;
        
		for (int ijet=0; ijet<ijetnum; ijet++){
            if(ijet_pt[ijet]>140 or ijet_pt[ijet]<0) continue;
            if(fabs(ijet_eta[ijet])>3) continue;
            leading_pt=0;
			for (int part=0; part<partnum; part++){
                
				dr = getr(ijet_eta[ijet],p_eta[part],ijet_phi[ijet],p_phi[part]);
				if ( dr > 0.7 ) continue;
				if( p_pt[part]>leading_pt) {
					leading_pt = p_pt[part];
					leading_index = part; 
				}
			}//part
            nijet++;

            tagged_id = fabs( p_id[leading_index] );
            if( tagged_id==1 or tagged_id==2 or tagged_id==3 ){
                nqtaggedijet++;
                
                qpt->Fill(ijet_pt[ijet], p_pt[leading_index]);
                qeta->Fill(ijet_eta[ijet], p_eta[leading_index]);
                qphi->Fill(ijet_phi[ijet], p_phi[leading_index]);
                qdr->Fill( getr(ijet_eta[ijet],p_eta[leading_index],ijet_phi[ijet],p_phi[leading_index]) );

                if( ijet_pt[ijet]>10 and ijet_pt[ijet]<40 and p_pt[leading_index]>20 and p_pt[leading_index]<40) areaijet[0]++;

            }
            else if( tagged_id==21 ){
                ngtaggedijet++;

                gpt->Fill(ijet_pt[ijet], p_pt[leading_index]);
                geta->Fill(ijet_eta[ijet], p_eta[leading_index]);
                gphi->Fill(ijet_phi[ijet], p_phi[leading_index]);
                gdr->Fill( getr(ijet_eta[ijet],p_eta[leading_index],ijet_phi[ijet],p_phi[leading_index]) );

                if( ijet_pt[ijet]>10 and ijet_pt[ijet]<40 and p_pt[leading_index]>20 and p_pt[leading_index]<40) areaijet[1]++;

            }
            else{
                netaggedijet++;

                ept->Fill(ijet_pt[ijet], p_pt[leading_index]);
                eeta->Fill(ijet_eta[ijet], p_eta[leading_index]);
                ephi->Fill(ijet_phi[ijet], p_phi[leading_index]);
                edr->Fill( getr(ijet_eta[ijet],p_eta[leading_index],ijet_phi[ijet],p_phi[leading_index]) );

                if( ijet_pt[ijet]>10 and ijet_pt[ijet]<40 and p_pt[leading_index]>20 and p_pt[leading_index]<40) areaijet[2]++;

            }

		}//ijet
    }//ien
    cout<<"all jets in condition: "<<nijet<<endl
        <<"light quark tagged: "<<nqtaggedijet<<"   jet pT 10~40, lp pT 20~40: "<<areaijet[0]<<endl
        <<"gluon quark tagged: "<<ngtaggedijet<<"   jet pT 10~40, lp pT 20~40: "<<areaijet[1]<<endl
        <<"else: "<<netaggedijet<<"   jet pT 10~40, lp pT 20~40: "<<areaijet[2]<<endl;
    
    TString outpath;
    if(eventmode == 1) outpath="./typesofjet2_fromgg2qqbar";
    else if( eventmode ==2) outpath = "./typesofjet2_fromgg2gg";
    else outpath="./typesofjet2_fromall";

    TString Rsize;
	Rsize.Form("R%02d",Rindex+2);
    outpath.Append(Rsize);
    outpath.Append(".root");
    TFile *outfile = new TFile(outpath, "RECREATE");
    
    qpt->Write(); qeta->Write(); qphi->Write(); qdr->Write();
    gpt->Write(); geta->Write(); gphi->Write(); gdr->Write();
    ept->Write(); eeta->Write(); ephi->Write(); edr->Write();

    outfile->Write();
    outfile->Close();
    return;
}
