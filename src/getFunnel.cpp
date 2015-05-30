#include "getFunnel.h"
void getFunnel(long line, Rcpp::IntegerVector& mother, Rcpp::IntegerVector& father, int* funnel, int nFounders)
{
	int currentLine = line;
	while(father(currentLine) == mother(currentLine))
	{
		currentLine = father(currentLine)-1;
	}
	int motherRow = mother(currentLine), fatherRow = father(currentLine);
	if(nFounders == 2)
	{
		funnel[0] = motherRow;
		funnel[1] = fatherRow;
		return;
	}
	int mmRow = mother(motherRow), mfRow = father(motherRow), fmRow = mother(fatherRow), ffRow = father(fatherRow);
	if(nFounders == 4)
	{
		funnel[0] = mmRow;
		funnel[1] = mfRow;
		funnel[2] = fmRow;
		funnel[3] = ffRow;
		return;
	}
	if(nFounders >= 8)
	{
		int mmmRow = mother(mmRow), mmfRow = father(mmRow), mfmRow = mother(mfRow), mffRow = father(mfRow);
		int fmmRow = mother(fmRow), fmfRow = father(fmRow), ffmRow = mother(ffRow), fffRow = father(ffRow);
		funnel[0] = mmmRow;
		funnel[1] = mmfRow;
		funnel[2] = mfmRow;
		funnel[3] = mffRow;
		funnel[4] = fmmRow;
		funnel[5] = fmfRow;
		funnel[6] = ffmRow;
		funnel[7] = fffRow;
		if(nFounders == 16)
		{
		  int mmmmRow = mother(mmmRow), mmmfRow=father(mmmRow), mmfmRow=mother(mmfRow), mmffRow=father(mmfRow);
		  int mfmmRow=mother(mfmRow), mfmfRow=father(mfmRow), mffmRow=mother(mffRow), mfffRow=father(mffRow);
		  int fmmmRow=mother(fmmRow), fmmfRow=father(fmmRow), fmfmRow=mother(fmfRow), fmffRow=father(fmfRow);
		  int ffmmRow=mother(ffmRow), ffmfRow=father(ffmRow), fffmRow=mother(fffRow), ffffRow=father(fffRow);
		  funnel[0] = mmmmRow;
		  funnel[1] = mmmfRow;
		  funnel[2] = mmfmRow;
		  funnel[3] = mmffRow;
		  funnel[4] = mfmmRow;
		  funnel[5] = mfmfRow;
		  funnel[6] = mffmRow;
		  funnel[7] = mfffRow;
		  funnel[8] = fmmmRow;
		  funnel[9] = fmmfRow;
		  funnel[10] = fmfmRow;
		  funnel[11] = fmffRow;
		  funnel[12] = ffmmRow;
		  funnel[13] = ffmfRow;
		  funnel[14] = fffmRow;
		  funnel[15] = ffffRow;
		}
	}
	else
	{
		throw std::runtime_error("Input nFounders of getFunnel must be 2, 4, 8 or 16");
	}
}