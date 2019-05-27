#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <ctime>
#include <random>
using namespace std;

//乱数とLambelt余弦則に必要な角度の生成クラス
class Rambelt {
private:  
    mt19937 mt();
public:
	double r;
	Rambelt() = default;
	Rambelt(double _r) : r(_r) {}

	double R() {
    	srand(time(0));
    	mt19937 mt(rand());
    	uniform_real_distribution<double> rand001(0.01, 1.00);   // [0.01 1.00] 範囲の一様乱数
		r = rand001(mt);
		return r;
	};
	double omega() {
    	srand(time(0));
    	mt19937 mt(rand());
    	uniform_real_distribution<double> rand001(0.01, 1.00); 
		r = rand001(mt);
		return 2 * 3.141592653589793 * r;
	};
	double phi() {
    	srand(time(0));
    	mt19937 mt(rand());
    	uniform_real_distribution<double> rand001(0.01, 1.00);
		r = rand001(mt);
		return acos(sqrt(1 - r));
	};
};

int main(void) {

	//パラメータ読み込み部分
	cout << "Start parameter reading..."  << endl;
	ifstream stream("data.csv");
	if (!stream) {
		cout << "can't read the file. Exit the program with Enter." << endl;
		cin.get();
		return 0;
	}

	string line;
	double data[23][11];
	const string delim = ",";

	int row = 0;
	int col;
	while (getline(stream, line)) {
		col = 0;
		// delimを区切り文字として切り分け、doubleに変換してdata[][]に格納
		for (string::size_type spos, epos = 0;
			(spos = line.find_first_not_of(delim, epos)) != string::npos;) {
			string token = line.substr(spos, (epos = line.find_first_of(delim, spos)) - spos);
			data[row][col++] = stod(token);
		}
		++row;
	}

	//読み込めたか確認
	for (row = 0; row < 23; ++row) {
		for (col = 0; col < 11; ++col) {
			cout << data[row][col] << " ";
		}
		cout << endl;
	}

	cout << "Reading complete. start MCmethod..." << endl;

	//モンテカルロ法開始
    clock_t start = clock();//計測開始
	int PN[23] = {};//壁面が吸収した粒子の数を数える変数
	int n[23]= {};//面番号
	int sep[23]= {};//分割数
	int dc[23] = {};//方向特性（xyz空間における面の向きを数字に対応させている）
	int f = 0;//反射してきた面情報変数
	int h = 0;//カウンタ変数
	double rc = 0;//壁面到達時の判定に使う変数
	double e[23] = {}; //放射率
	double e2[23] = {};	//放射率+指向性反射率

	//判定に必要なパラメータの事前代入
	for (h = 0; h < 23; ++h) {
		n[h] = (int)data[h][0]; 
		sep[h] = (int)data[h][1];
		e[h] = data[h][2];
		e2[h] = data[h][3];
		dc[h] = (int)data[h][4];
	}

	//計算開始
	for (h = 0; h < 23; ++h) {

		//データ入力、開始地点(x.y.z)
		double x0 = data[h][5];
		double y0 = data[h][6];
		double z0 = data[h][7];
		//面の大きさ(したがってどれかは0になるはず)
		double lx = data[h][8];
		double ly = data[h][9];
		double lz = data[h][10];
		//反射してきた面情報変数初期化
		f = dc[h];
		cout << "No."<< n[h] << " parameter reading..." << endl;
		cout << x0 << " "<< y0 << " "<< z0 << " "<< lx << " " << ly << " " << lz << endl;
		cout <<"Start caliculating..." << endl;
		for (long i = 1; i < sep[h]; i++) {
			for (long j = 1; i < sep[h]; i++) {
				for (long k = 1; i < sep[h]; i++) {
					//初期位置代入
					double x = x0 + lx * i / (sep[h] + 1);
					double y = y0 + ly * j / (sep[h] + 1);
					double z = z0 + lz * k / (sep[h] + 1);
					//乱数による速度、方向生成
					Rambelt test;
					double Vx = sin(test.phi()) * cos(test.omega());
					double Vy = sin(test.phi()) * sin(test.omega());
					double Vz = cos(test.phi());

					//速度方向調整(体系内向かうように方向を変える)
					if (dc[f] == 1) {
						Vx = Vz;
						Vy = Vx;
						Vz = Vy;
						}
					else if(dc[f] == -1) {
						Vx = -1 * Vz;
						Vy = Vx;
						Vz = Vy;
						}
					else if(dc[f] == 2) {
						Vy = Vz;
						Vx = Vx;
						Vz = Vy;
						}
					else if(dc[f] == -2) {
						Vy = -1 * Vz;
						Vx = Vx;
						Vz = Vy;
						}
					else if(dc[f] == 3) {
						Vz = Vz;
						Vx = Vx;
						Vy = Vy;
						}
					else {
						Vz = -1 * Vz;
						Vx = Vx;
						Vy = Vy;
						}

					//粒子移動(壁面で反射or乱反射した場合はここに戻る)
					moving:
					x = x + Vx;
					y = y + Vy;
					z = z + Vz;

					//壁到達判定
					if(x < 0) {
						x = 0;

						if(y > 340) {
							f = 23;
							}
						else if(y < 20) {
							f = 24;
							}
						else if(y > 20 && y < 340 && z < 60) {
							f = 25;
							}
						else if(y > 20 && y < 340 && z > 440) {
							f = 26;
							}
						else {
							f = 5;
							}
						}
					else if(z > 501 && x < 75) {
						x = 75;

						if(z > 900) {
							f = 18;
							}
						else {
							f = 14;
							}
						}
					else if(x > 150) {
						x = 150;

						if(z < 500) {
							f = 11;
							}
						else if(z > 900) {
							f = 20;
							}
						else {
							f = 16;
							}
						}
					else if(y < 0) {
						y = 0;

						if(z < 500) {
							if(x < 75) {
								f = 7;
								}
							else {
								f = 12;
								}
							}
						else if(z > 900) {
							f = 21;
							}	
						else {
							f = 17;
							}
						}
					else if(y > 360) {
						y = 360;

						if(z < 500) {
							if(x < 75) {
								f = 6;
								}
							else {
								f = 10;
								}
							}
						else if(z > 900) {
							f = 19;
							}
						else {
							f = 15;
							}
						}
					else if(z < 0) {
						z = 0;

						if(x < 75) {
							f = 8;
							}
						else {
							f = 13;
							}
						}	
					else if(z > 500 && x < 75) {
						z = 500;
						f = 9;
						}
					else if(z > 1300) {
						z = 1300;
						f = 22;
						}
					else {
						goto moving;//禁断のgoto文、movingに戻る
						}

					//壁面到達時の挙動操作
					rc = test.R();//壁面衝突時の挙動判定の為の乱数生成			
     				//吸収する場合
					if(rc < e[f]) { 
							PN[f] += 1;//吸収した粒子の数+1
							cout << "particle is absorbed."<< endl;
						}
					//乱反射する場合
					else if(rc > e2[f]) {
						Vx = sin(test.phi()) * cos(test.omega());
						Vy = sin(test.phi()) * sin(test.omega());
						Vz = cos(test.phi());
						if (dc[f] == 1) {
							Vx = Vz;
							Vy = Vx;
							Vz = Vy;
							}
						else if(dc[f] == -1) {
							Vx = -1 * Vz;
							Vy = Vx;
							Vz = Vy;
							}
						else if(dc[f] == 2) {
							Vy = Vz;
							Vx = Vx;
							Vz = Vy;
							}
						else if(dc[f] == -2) {
							Vy = -1 * Vz;
							Vx = Vx;
							Vz = Vy;
							}
						else if(dc[f] == 3) {
							Vz = Vz;
							Vx = Vx;
							Vy = Vy;
							}
						else if(dc[f] == -3) {
							Vz = -1 * Vz;
							Vx = Vx;
							Vy = Vy;
							}
						goto moving;//吸収時以外はmoving（粒子の移動）に戻る	
						}

					//指向性反射する場合
					else {                         
						if(dc[f] == 1 || dc[f] == -1) {
							Vx = -1 * Vx;
							}
						else if(dc[f] == 2 || dc[f] == -2) {
							Vy = -1 * Vy;
							}
						else {
							Vz = -1 * Vz;
							}
						goto moving;//吸収時以外はmoving（粒子の移動）に戻る			
						}
				}
			}
		}
		cout << "No." << n[h] << " finished. NextFace..." << endl;
	}
	clock_t end = clock();	//計測終了
	//粒子総数確認
	int p = 0;//粒子総数を数える変数
	for( h = 0; h < 23; ++h){
			p += PN[h];
		}	
	cout <<"Total number of particles = " << PN[h] << endl;//この値がsep(分割数)の合計と等しくなっていれば少なくともエネルギーは保存されている

    //計算結果出力（単位:個）
    ofstream ofs("result.csv");
    if (!ofs){
        cout << "can't read the file. Exit the program with Enter." << endl;
        cin.get();
        return 0;
    	}

	for (h = 0; h < 23; ++h) {
    		ofs << n[h] << "," << PN[h] << "," << endl;
		}
	ofs.close();
	cout << "Wrinting complete." << endl;
    //プログラム終了と計算時間表示
	cout << "Caliculation time is " << (double)(end - start) / CLOCKS_PER_SEC << " seconds.\nExit the program with Enter.\n" << endl;
	cin.get();
	return 0;
}
