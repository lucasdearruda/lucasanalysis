std::vector<int> runsForward = {371, 373, 380, 380};


void test(){
	for (Int_t intv = 0; intv + 1 < runsForward.size(); intv += 2)
	{
		int startRun = runsForward[intv];
		int endRun = runsForward[intv + 1];
		for (int i = startRun; i <= endRun; ++i)
		{
			cout<<i<<endl;
		}
	}
	return;
}
