/***********************************************************************
*  This code is part of PolyVest.
*
*  Copyright (C) 2013, 2016 Cunjing Ge, Institute of Software, Chinese 
*  Academy of Sciences, Beijing, China. All rights reserved. 
*  E-mail: <gecj@ios.ac.cn>.
*
*  PolyVest is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
* 
*  PolyVest is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with PolyVest. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#include <fstream>
#include "vol.h"

using namespace vol;

ifstream file;

double get_num(){
	char ch = ' ';
	double tmp_num = 0;
	bool isnumber = false;
	bool tmp_neg = false;
	bool decimal = false;
	int dec_len = 0;
	for (; (!file.eof()) && (isspace(ch)); file.get(ch));
	while (!file.eof()){
		if (isdigit(ch)){
			if (decimal) dec_len++;
			tmp_num = tmp_num * 10 + ch - 48;
			isnumber = true;
		}else if (ch == '-') tmp_neg = true;
		else if (ch == '.')	decimal = true;
		else if (isspace(ch)) break;
		else{
			cout << "error: unknown identifier \'" << ch << "\'." << endl;
			exit(1);
		}
		file.get(ch);
	}
	if (!isnumber){
		cout << "error: input file is incomplete, please check the size of matrix." << endl;
		exit(1);
	}
	if (tmp_neg) tmp_num = -tmp_num;
	if (dec_len > 0) tmp_num /= pow((double)10, dec_len);
	return tmp_num;
}

int main(int argc, char *argv[]){

	cout << endl << "--------------- PolyVest ----------------" << endl;
	cout << "If you have any questions or if you have found some bugs," << endl << "please contact me at <gecj@ios.ac.cn>." << endl;
	cout << endl << "=========================================" << endl;

	if (argc <= 2 || argc >= 5){
		cout << "error: invalid arguments." << endl;
		cout << "USAGE: " << argv[0] << " <input-file>  <step-size-coef> [output-file]" << endl;
		return 1;
	}

	file.open(argv[1]);

	if (!file.is_open()){ cout << "Cant open file." << endl; return 1; }

	int rows, cols;
	rows = get_num();
	cols = get_num();
	
	polytope p(rows, cols);

	cout << "Hyperplanes: " << rows << endl << "Dimensions:" << cols << endl;

	for (int i = 0; i < rows; i++){
		double t = get_num();
		p.vecb(t, i);
		for (int j = 0; j < cols; j++){
			t = get_num();
			p.matA(t, i, j);
		}
	}

	//p.msg_off = true;
	p.check_planes_off = true;
	cout << endl << "============= Preprocessing =============" << endl;
	p.Preprocess();

	cout << endl << "=============== Sampling ================" << endl;
	cout << endl << "Volume of polytope (estimate): " << p.EstimateVol(atoi(argv[2])) << endl;
	ofstream outfile;
	if (argc == 3){
		outfile.open("PolyVest.result", ios::app);
	}else{
		outfile.open(argv[3]);
	}
	outfile << p.Volume() << endl;
	outfile.close();

/*
	ofstream outfile(argv[3]);
	for (int i = 0; i < ((argc == 5) ? atoi(argv[4]) : 100); i++){
		if (i % 20 == 0) cout << i << endl;
		p.EstimateVol(atoi(argv[2]));
		outfile << p.Volume() << endl;
	}
	outfile.close();
*/
	//print matrix A and vector b which applied affine transformation
/*
	for (int i = 0; i < rows; i++){
		cout << p.vecb(i) << ' ';
		for (int j = 0; j < cols; j++) cout << p.matA(i, j) << ' ';
		cout << endl;
	}
*/

	return 0;
}
