#include "disjoint_set.h"

set* set::find(void) {									// Declares a function 'find' from the set class and returns a pointer 

	if (parent != this) parent = parent->find();		// Checks if 
	return parent;
}

void set::unite(set* other) {

	set* thisroot = this->find();
	set* otherroot = other->find();

	if (otherroot == thisroot) return;

	int thisrank = thisroot->get_rank();
	int otherrank = otherroot->get_rank();

	if (thisrank < otherrank) {

		thisroot->set_parent(otherroot);
	}
	else if (thisrank > otherrank) {

		otherroot->set_parent(thisroot);
	}
	else {

		otherroot->set_parent(thisroot);
		thisroot->upgrade_rank();
	}
}
