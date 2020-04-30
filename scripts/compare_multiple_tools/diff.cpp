#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream> 

using namespace std;

void separate_by_comma(string s, vector <string>& v) {
	v.clear();
	if (s.length() == 0)
		return;

	std::stringstream ss(s);
	
	while (ss.good()) {
		string substr;
		getline(ss, substr, ',');
		v.push_back(substr);
	}
}

// b: before
// a: after
void get_difference(string b, string a) {
	vector <string> vb;
	vector <string> va;

	separate_by_comma(b, vb);
	separate_by_comma(a, va);

	sort(vb.begin(), vb.end());
	sort(va.begin(), va.end());

	vector <string> b_a(vb.size());
	vector <string> a_b(va.size());

	auto it_b_a = std::set_difference(vb.begin(), vb.end(), va.begin(), va.end(), b_a.begin());
	auto it_a_b = std::set_difference(va.begin(), va.end(), vb.begin(), vb.end(), a_b.begin());

	b_a.resize(it_b_a - b_a.begin());
	a_b.resize(it_a_b - a_b.begin());

	cout << "\t";
	if (b_a.begin() == b_a.end())
		cout << "-";
	for (auto it = b_a.begin(); it != b_a.end(); ++it)
		cout << *it << ",";

	cout << "\t";
	if (a_b.begin() == a_b.end())
		cout << "-";
	for (auto it = a_b.begin(); it != a_b.end(); ++it)
		cout << *it << ",";
}

int main(int argc, char** argv) {
	--argc;

	if (argc < 2)
		cerr << "Usage: " << argv[0] << " FILE1 FILE2 [FILE3 FILE4]\n";

    ifstream fin1;
    ifstream fin2;
    ifstream fin3;
    ifstream fin4;

    fin1.open(argv[1]);
    fin2.open(argv[2]);
    fin3.open(argv[3]);
	fin4.open(argv[4]);

    string a, b, c, rs;
    int d;
    map <string, vector<int> > mymap;
    map <string, pair<string, string> > rsmap;
    while (fin1 >> a >> b >> c >> d >> rs) {
        mymap[a+'-'+b+'-'+c] = {d, 0, 0, 0};
		rsmap[a+'-'+b+'-'+c] = {rs, ""};
    }

    while (fin2 >> a >> b >> c >> d >> rs) {
        if (mymap.find(a+'-'+b+'-'+c) != mymap.end()) {
            mymap[a+'-'+b+'-'+c][1] = d;
			rsmap[a+'-'+b+'-'+c].second = rs;
        }
        else {
            mymap[a+'-'+b+'-'+c] = {0, d, 0, 0};
			rsmap[a+'-'+b+'-'+c] = {"", rs};
        }
    }

    while (fin3 >> a >> b >> c >> d) {
        if (mymap.find(a+'-'+b+'-'+c) != mymap.end()) {
            mymap[a+'-'+b+'-'+c][2] = d;
        }
    }

	int ce2b;
    while (fin4 >> a >> ce2b >> c >> d) {
		b = to_string(ce2b+1);
        if (mymap.find(a+'-'+b+'-'+c) != mymap.end()) {
            mymap[a+'-'+b+'-'+c][3] = d;
        }
    }

    map <string, pair<int, int> >::iterator it;
    int increased = 0, decreased = 0, not_changed = 0, disapeared = 0, newly_added = 0;

	auto rsit = rsmap.begin();
    for (auto it = mymap.begin(); it != mymap.end(); ++it) {
        int bb = it->second[0];
        int aa = it->second[1];

		int ciri = it->second[2];
		int ce2 = it->second[3];

        if (bb - aa == 0) {
            not_changed++;
            cout << "[NC]\t";
        }
        else if (bb < aa) {
            if (bb == 0) {
                newly_added++;
                cout << "[NA]\t";
            }
            else {
                increased++;
                cout << "[++]\t";
            }
        }
        else {
            if (aa == 0) {
                disapeared++;
                cout << "[DI]\t";
            }
            else {
                decreased++;
                cout << "[--]\t";
            }
        }

        cout << it->first << "\t" << bb << "\t" << aa;
		if (argc > 2)
        	cout << "\t" << ciri << "\t" << ce2;

		get_difference(rsit->second.first, rsit->second.second);

		cout << endl;

		++rsit;

    }

	fin1.close();
	fin2.close();
	fin3.close();
	fin4.close();

    return 0;
}
