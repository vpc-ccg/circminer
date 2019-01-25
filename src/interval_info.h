#ifndef __INTERVAL_INFO_H__
#define __INTERVAL_INFO_H__

#include <stdint.h>
#include <vector>

using namespace std;

template <class T> class IntervalInfo {
public:
	uint32_t spos;
	uint32_t epos;
	uint32_t max_end;
	uint32_t min_end;
	uint32_t max_next_exon;

	vector<T> seg_list;

public:
	IntervalInfo(const T& seg) {
		spos = seg.start;
		epos = seg.end;
		seg_list.clear();
		seg_list.push_back(seg);
	}
	
	IntervalInfo(uint32_t s, uint32_t e, const vector<T>& seg_l) {
		spos = s;
		epos = e;
		seg_list.clear();
		seg_list.resize(seg_l.size());
		for (int i = 0; i < seg_l.size(); i++) {
			seg_list[i] = seg_l[i];
		}
	}

	IntervalInfo(uint32_t s, uint32_t e, const vector<T>& seg_l, const T& seg) : IntervalInfo(s, e, seg_l) {
		seg_list.push_back(seg);
	}

	IntervalInfo(const IntervalInfo& ii) {
		spos = ii.spos;
		epos = ii.epos;
		seg_list.clear();
		seg_list.resize(ii.seg_list.size());
		for (int i = 0; i < ii.seg_list.size(); i++) {
			seg_list[i] = ii.seg_list[i];
		}
	}

	IntervalInfo& operator = (const IntervalInfo& ii) {
		if (this == &ii)	// self assignment prevention
			return *this;

		spos = ii.spos;
		epos = ii.epos;
		seg_list.clear();
		seg_list.resize(ii.seg_list.size());
		for (int i = 0; i < ii.seg_list.size(); i++) {
			seg_list[i] = ii.seg_list[i];
		}

		return *this;
	}
	
	void push_back(const T& seg) {
		seg_list.push_back(seg);
	}
};

#endif	//__INTERVAL_INFO_H__
