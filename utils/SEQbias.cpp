#include "SEQbias.hpp"


void ReadSeqBias(string gzip_file, int32_t& contextLeft, int32_t& contextRight, vector<int32_t>& orders, 
	vector<int32_t>& shifts, vector<int32_t>& widths, Eigen::MatrixXd& Probs)
{
	// file handle
	boost::iostreams::filtering_istream fpin;
	fpin.push(boost::iostreams::gzip_decompressor());
	fpin.push(boost::iostreams::file_source(gzip_file, std::ios_base::in | std::ios_base::binary));

	// read file content
	int32_t contextLength;
	fpin.read((char*)&contextLength, sizeof(int32_t));
	fpin.read((char*)&contextLeft, sizeof(int32_t));
	fpin.read((char*)&contextRight, sizeof(int32_t));
	assert(contextLength == contextLeft + contextRight + 1);

	orders.resize(contextLength);
	shifts.resize(contextLength);
	widths.resize(contextLength);
	fpin.read((char*)&orders[0], contextLength*sizeof(int32_t));
	fpin.read((char*)&shifts[0], contextLength*sizeof(int32_t));
	fpin.read((char*)&widths[0], contextLength*sizeof(int32_t));

	typename Eigen::MatrixXd::Index prows, pcols;
	fpin.read((char*)&prows, sizeof(typename Eigen::MatrixXd::Index));
	fpin.read((char*)&pcols, sizeof(typename Eigen::MatrixXd::Index));

	Probs=Eigen::MatrixXd::Zero(prows, pcols);
	fpin.read((char*)(Probs.data()), prows*pcols*sizeof(typename Eigen::MatrixXd::Scalar));
};


SeqBiasModel_t::SeqBiasModel_t(string obs5_gzip_file, string obs3_gzip_file, string exp5_gzip_file, string exp3_gzip_file)
{
	time_t CurrentTime;
	string CurrentTimeStr;
	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout<<"["<<CurrentTimeStr.substr(0, CurrentTimeStr.size()-1)<<"] "<<"Loading sequence bias model."<<endl;

	// reading obs5
	ReadSeqBias(obs5_gzip_file, contextLeft, contextRight, orders, shifts, widths, Obs5Probs);

	// reading obs3
	int32_t tmp_contextLeft;
	int32_t tmp_contextRight;
	vector<int32_t> tmp_orders, tmp_shifts, tmp_widths;
	ReadSeqBias(obs3_gzip_file, tmp_contextLeft, tmp_contextRight, tmp_orders, tmp_shifts, tmp_widths, Obs3Probs);
	assert(tmp_contextLeft == contextLeft && tmp_contextRight == contextRight);
	assert(tmp_orders.size() == orders.size() && tmp_widths.size() == widths.size() && tmp_shifts.size() == shifts.size());
	for (int32_t i = 0; i < contextLeft + contextRight + 1; i++)
		assert(tmp_orders[i] == orders[i] && tmp_widths[i] == widths[i] && tmp_shifts[i] == shifts[i]);

	// reading exp5
	ReadSeqBias(exp5_gzip_file, tmp_contextLeft, tmp_contextRight, tmp_orders, tmp_shifts, tmp_widths, Exp5Probs);
	assert(tmp_contextLeft == contextLeft && tmp_contextRight == contextRight);
	assert(tmp_orders.size() == orders.size() && tmp_widths.size() == widths.size() && tmp_shifts.size() == shifts.size());
	for (int32_t i = 0; i < contextLeft + contextRight + 1; i++)
		assert(tmp_orders[i] == orders[i] && tmp_widths[i] == widths[i] && tmp_shifts[i] == shifts[i]);

	// reading exp3
	ReadSeqBias(exp3_gzip_file, tmp_contextLeft, tmp_contextRight, tmp_orders, tmp_shifts, tmp_widths, Exp3Probs);
	assert(tmp_contextLeft == contextLeft && tmp_contextRight == contextRight);
	assert(tmp_orders.size() == orders.size() && tmp_widths.size() == widths.size() && tmp_shifts.size() == shifts.size());
	for (int32_t i = 0; i < contextLeft + contextRight + 1; i++)
		assert(tmp_orders[i] == orders[i] && tmp_widths[i] == widths[i] && tmp_shifts[i] == shifts[i]);

	time(&CurrentTime);
	CurrentTimeStr=ctime(&CurrentTime);
	cout << "[" << CurrentTimeStr.substr(0, CurrentTimeStr.size()-1) << "] " << "Finish.\n";
};