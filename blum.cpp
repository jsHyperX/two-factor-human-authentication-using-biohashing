#include<bits/stdc++.h>
#define _ ios_base::sync_with_stdio(0);cin.tie(0);
using namespace std;
const int N=100001;
typedef unsigned long long ll;
int prime[N];                   // prime[i] tells if number i is prime or not
vector<ll> pt;                  // this vector stores a list of primes which is later used to pick up two blum primes
ll gcd(ll a,ll b) {
    if(!b) return a;
    return gcd(b,a%b);
}
ll fast_exp(ll base,ll exp,ll m) {      // fast exponentiation in O(logn)
    ll res=1;
    while(exp>0){
        if(exp&1)  res = (res * base)%m;
        base = (base * base)%m;
        exp= exp >> 1;
    }
    return res;
}
/*****************************************************************************
    Carmichael function for an integer is given by lambda(p*q) = lcm(p-1,q-1);
    so, lambda(p*q) = ((p-1)*(q-1))/gcd(p-1,q-1)
****************************************************************************/
ll carmichael(ll a,ll b) {
    a--;b--;
    ll num=((1LL)*a)*b;
    ll den=gcd(a,b);
    return num/den;
}
/*********************************************************************************
		generating prime numbers using sieve of Eratosthenes
**********************************************************************************/
void gen_p() {
    prime[0]=prime[1]=0;
    for(int i=2;i<N;i++) {
        if(prime[i]) {
            pt.push_back(i);            //storing the prime number i in vector pt
            for(int j=i+i;j<N;j+=i) prime[j]=0;
        }
    }
    //for(int i=0;i<pt.size();i++) cout << pt[i] << " ";
    //cout << endl;
}
/***********************************************************************************
        the blum blum PRBG involves the following recurrence x_(n+1)=(x_n)^2 mod (m)
        where m is a prime such that m=p*q,p and q are primes such that p=3(mod 4),
        q=3(mod 4) and gcd(phi(p-1),phi(q-1)) is small
***********************************************************************************/
void gen(ll seed,ll i,ll p,ll q) {
    ll ans = fast_exp(seed,fast_exp(2,i,carmichael(p,q)),(1LL)*p*q);
    cout << ((ans&1)?1:0) << " ";
    //cout << ans << " ";
}
void solve() {
    ll n,m,p[2],s;
	cin >> m;
	memset(prime,1,sizeof(prime));              //initializing the array prime with all 1
	cout << "generating prime numbers......" << endl;
	gen_p();            // generating primes
    ll si=pt.size()-1;
    cout << "................................." << endl;
    int count=0;
    cout << "generating two Blum prime numbers......" << endl;
    srand(time(NULL));
    while(count<2) {
        int a=rand()%si;                // generating a random index in the range of size of array pt
        ll temp=pt[si-a];               // picking up a prime number using the generated index
        if(temp%4==3) { p[count]=temp;count++; }            // checking if the prime obtained is blum or not
    }
    //cout << p[0] << " " << p[1] << " ";
    cout << "blum primes are generated.............." << endl;
    n = ((1LL)*p[0])*p[1];
    cout << "................................" << endl;
    cout << "................................" << endl;
    cout << "generating seed............." << endl;
    s=rand()%n;                         // generating a random seed within n
    cout << "seed generated......." << endl;
    cout << "now generating the PR Bit Sequence.............." << endl;
    cout << "here you go:" << " ";
    for(int i=2;i<=m;i++)  gen(s,i,p[0],p[1]);
}
int main() {
    int t;cin >> t;
    while(t--) {
        solve();
        cout << "\n";
    }
	return 0;
}

