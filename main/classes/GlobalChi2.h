//GLoal Chi2 operator for simutaneous fit
struct GlobalChi2 {
   const  ROOT::Math::IMultiGenFunction* fChi2_1;
   const  ROOT::Math::IMultiGenFunction* fChi2_2;


	int iparPass[12] = {0,
						1,
						2,
						3,
						4,
						5,
						6,
						7,
						8,
						9,
						10,
						11,
					};

	int iparPassFail[24] = {0,
						1,
						2,
						3,
						4,
						5,
						6,
						7,	
						8,
						9,
						10,
						11,
						12,
						13,
						14,
						15,
						16,
						17,
						18,
						19,
						20,
						21,
						22,
						23,
					};

   double operator() (const double* par) const {
      double p1[12];
      for (int i = 0; i < 12; ++i) p1[i] = par[iparPass[i]];

      double p2[24];
      for (int i = 0; i < 24; ++i) p2[i] = par[iparPassFail[i]];

      return (*fChi2_1)(p1) + (*fChi2_2)(p2);
   }

   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2) :
      fChi2_1(&f1), fChi2_2(&f2) {}
};