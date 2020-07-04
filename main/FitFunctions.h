
//-------------------------------------
// Fit functions for Invariant Mass
//-------------------------------------

class FitFunctions
{
public:
	class Primary
	{
	public:
		//Gaussian function
		static Double_t Gaus(Double_t *x, Double_t *par)
		{
			//par[0] = height
			//par[1] = position
			//par[2] = sigma
			Double_t gaus = par[0]*TMath::Gaus(x[0],par[1],par[2],1);
			return gaus;
		}

		//Polynomial function
		static Double_t Pol1(Double_t *x, Double_t *par)
		{
			//par[0] = b
			//par[1] = a
			Double_t pol = par[0] + par[1]*x[0];
			return pol;
		}

		//Exponential function
		static Double_t Exp(Double_t *x, Double_t *par)
		{
			//par[0] = height
			//par[1] = width
			Double_t exp = par[0] * TMath::Exp(par[1]*x[0]);
			return exp;
		}

		//crystall ball function
		static Double_t CrystalBall(Double_t *x,Double_t *par)
		{
			//par[0] = alpha
			//par[1] = n
			//par[2] = mean
			//par[3] = sigma
			//par[4] = Yield
			Double_t t = (x[0]-par[2])/par[3];
			if (par[0] < 0) t = -t;
			Double_t absAlpha = fabs((Double_t)par[0]);
			if (t >= -absAlpha)
			{
				return par[4]*exp(-0.5*t*t);
			}
			else
			{
				Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
				Double_t b = par[1]/absAlpha - absAlpha;
				return par[4]*(a/TMath::Power(b - t, par[1]));
			}
		}
	};
	class Merged
	{
	public:

		//OLD

		//Fit function for signal for Invariant Mass
		static Double_t Signal_InvariantMassAll(Double_t *x, Double_t *par) {
			return FitFunctions::Primary::Gaus(x,par) + FitFunctions::Primary::CrystalBall(x, &par[3]);
		}

		//Fit function for background for Invariant Mass
		static Double_t Background_InvariantMassAll(Double_t *x, Double_t *par) {
			return FitFunctions::Primary::Exp(x,par) + FitFunctions::Primary::Exp(x, &par[2]);
		}

		//Fit function for signal & background for Invariant Mass
		static Double_t FFit_InvariantMassAll(Double_t *x, Double_t *par) {
			return FitFunctions::Merged::Signal_InvariantMassAll(x,par) + FitFunctions::Merged::Background_InvariantMassAll(x, &par[8]);
		}

		//NEW

		//Fit function for signal & background for Invariant Mass
		static Double_t Pass_InvariantMass(Double_t *x, Double_t *par) {
			return FitFunctions::Merged::Signal_InvariantMassAll(x,par) + FitFunctions::Merged::Background_InvariantMassAll(x, &par[8]);
		}

		//Fit function for signal & background for Invariant Mass
		static Double_t Fail_InvariantMass(Double_t *x, Double_t *par) {
			return FitFunctions::Merged::Signal_InvariantMassAll(x,par) + FitFunctions::Merged::Background_InvariantMassAll(x, &par[8]);
		}

		//Fit function for signal & background for Invariant Mass
		static Double_t Both_InvariantMass(Double_t *x, Double_t *par) {
			return FitFunctions::Merged::Pass_InvariantMass(x,par) + FitFunctions::Merged::Fail_InvariantMass(x, &par[12]);
		}
	};
};