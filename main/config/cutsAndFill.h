//This files holds some functions used in macro.cpp for particle selection

//Return if is a accepted particle or no
bool applyCuts(double** quantities, int** types)
{
	//Assign variables for easy visualization
	double &ProbeMuon_Pt            = *quantities[0];
	double &ProbeMuon_Eta           = *quantities[1];
	double &ProbeMuon_Phi           = *quantities[2];
	double &TagMuon_Pt              = *quantities[3];
	double &TagMuon_Eta             = *quantities[4];
	double &TagMuon_Phi             = *quantities[5];
	double &InvariantMass           = *quantities[6];
	int &PassingProbeTrackingMuon   = *types[0];
	int &PassingProbeStandAloneMuon = *types[1];
	int &PassingProbeGlobalMuon     = *types[2];

	//Apply cuts
	if (TagMuon_Pt >= 7.0 && abs(TagMuon_Eta) <= 2.4)
		return true;

	return false;
}

void fillMassHistograms(Particle *MotherObj, double** quantities, int** types)
{
	//Assign variables for easy visualization
	double &ProbeMuon_Pt            = *quantities[0];
	double &ProbeMuon_Eta           = *quantities[1];
	double &ProbeMuon_Phi           = *quantities[2];
	double &TagMuon_Pt              = *quantities[3];
	double &TagMuon_Eta             = *quantities[4];
	double &TagMuon_Phi             = *quantities[5];
	double &InvariantMass           = *quantities[6];
	int &PassingProbeTrackingMuon   = *types[0];
	int &PassingProbeStandAloneMuon = *types[1];
	int &PassingProbeGlobalMuon     = *types[2];
	Particle &Muon = *MotherObj;
	
	//Filling histograms are below

	if (PassingProbeTrackingMuon)
		Muon.Tracker.hMass->Fill(InvariantMass);

	if (PassingProbeStandAloneMuon)
		Muon.Standalone.hMass->Fill(InvariantMass);

	if (PassingProbeGlobalMuon)
		Muon.Global.hMass->Fill(InvariantMass);

	Muon.All.hMass->Fill(InvariantMass);
}

//Fill histogram
void fillQuantitiesHistograms(Particle *MotherObj, double** quantities, int** types)
{
	//Assign variables for easy visualization
	double &ProbeMuon_Pt            = *quantities[0];
	double &ProbeMuon_Eta           = *quantities[1];
	double &ProbeMuon_Phi           = *quantities[2];
	double &TagMuon_Pt              = *quantities[3];
	double &TagMuon_Eta             = *quantities[4];
	double &TagMuon_Phi             = *quantities[5];
	double &InvariantMass           = *quantities[6];
	int &PassingProbeTrackingMuon   = *types[0];
	int &PassingProbeStandAloneMuon = *types[1];
	int &PassingProbeGlobalMuon     = *types[2];
	Particle &Muon = *MotherObj;

	//Filling histograms are below

	if (PassingProbeTrackingMuon)
		Muon.Tracker.fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);

	if (PassingProbeStandAloneMuon)
		Muon.Standalone.fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);

	if (PassingProbeGlobalMuon)
		Muon.Global.fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);

	Muon.All.fillHistograms(InvariantMass, TagMuon_Pt, TagMuon_Eta, TagMuon_Phi, ProbeMuon_Pt, ProbeMuon_Eta, ProbeMuon_Phi);
}