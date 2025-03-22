#ifndef DISJOINT_SET_H
#define DISJOINT_SET_H

class set {
private:

	set* parent;
	int rank;
public:

	set(){ rank = 0; parent = this; }
	~set(){}

	int get_rank() { return rank; }
	void set_parent (set* value) { parent = value; }
	void upgrade_rank () { ++rank; }

	set* find();
	void unite(set*);
};

#endif
