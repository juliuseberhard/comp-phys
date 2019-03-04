// &x means "memory address of variable x"
// *px means "value at memory address px"

#include <iostream>

using namespace std;

void f(int *x) {
	*x = *x + 1;
}

int main() {
  int x;
	x = 2;
	cout << "x = " << x << endl;  // x = 2
	f(&x);  // changes value of x within the memory
	        // without return statement in f
  cout << "x = " << x << endl;  // x = 3
	cout << "&x = " << &x << endl;  // prints the address of x

  int *py;  // empty address py, was not specified before
	          // instead, a default address is assigned:
	cout << py << endl;  // prints the default address
	
	return 0;
}
