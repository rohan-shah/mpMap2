template<int nFounders, bool infiniteSelfing> struct compressedProbabilities;
template<> struct compressedProbabilities<2, true>
{
	static const int nDifferentProbs = 2;
};
template<> struct compressedProbabilities<2, false>
{
	static const int nDifferentProbs = 5;
};
template<> struct compressedProbabilities<4, true>
{
	static const int nDifferentProbs = 3;
};
template<> struct compressedProbabilities<4, false>
{
	static const int nDifferentProbs = 18;
};
template<> struct compressedProbabilities<8, true>
{
	static const int nDifferentProbs = 3;
};
template<> struct compressedProbabilities<8, false>
{
	static const int nDifferentProbs = 46;
};
template<> struct compressedProbabilities<16, true>
{
	static const int nDifferentProbs = 4;
};
template<> struct compressedProbabilities<16, false>
{
	static const int nDifferentProbs = 95;
};
