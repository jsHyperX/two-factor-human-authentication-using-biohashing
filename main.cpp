#include<bits/stdc++.h>
#define _ ios_base::sync_with_stdio(0);cin.tie(0);
using namespace std;
typedef vector<float> vi;
typedef vector<vi> vii;
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
}
/***********************************************************************************
        the blum blum PRBG involves the following recurrence x_(n+1)=(x_n)^2 mod (m)
        where m is a prime such that m=p*q,p and q are primes such that p=3(mod 4),
        q=3(mod 4) and gcd(phi(p-1),phi(q-1)) is small
***********************************************************************************/
ll gen(ll seed,ll i,ll p,ll q) {
    ll ans = fast_exp(seed,fast_exp(2,i,carmichael(p,q)),(1LL)*p*q);
    return (ans&1)?1:0;
}
vi solve(int m) {
    vi ans;
    ll n,p[2],s;
	memset(prime,1,sizeof(prime));              //initializing the array prime with all 1
	gen_p();            // generating primes
    ll si=pt.size()-1;
    int count=0;
    srand(time(NULL));
    while(count<2) {
        int a=rand()%si;                // generating a random index in the range of size of array pt
        ll temp=pt[si-a];               // picking up a prime number using the generated index
        if(temp%4==3) { p[count]=temp;count++; }            // checking if the prime obtained is blum or not
    }
    n = ((1LL)*p[0])*p[1];
    s=rand()%n;
    for(int i=2;i<=m;i++) {
        ans.push_back(gen(s,i,p[0],p[1]));
    }
    return ans;
}
vi transform(vi a,vi b,vi c) {
    vi vec;
    float coeff;
    float num=0,denum=0;
    for(int i=0;i<a.size();i++) {
        num+=a[i]*b[i];denum+=b[i]*b[i];
    }
    coeff=num/denum;
    float mag=0;
    for(int i=0;i<a.size();i++) {
        float val=c[i]-(coeff*b[i]);
        mag+=val*val;
        vec.push_back(val);
    }
    mag=sqrt(mag);
    for(int i=0;i<vec.size();i++) vec[i]/=mag;
    return vec;
}
vii orthonormal(vii vec) {
    vii u;
    int rows=vec.size();
    u.push_back(vec[0]);
    for(int i=1;i<rows;i++) {
        vi v = vec[i];
        int k=i;
        while(k) {
            v=transform(vec[i],u[k-1],v);k--;
        }
        u.push_back(v);
    }
    return u;
}
/**************************************************************************************************
                function to take facial data as a input in the form of a 2-D matrix
**************************************************************************************************/
void inp(vii matrix,int n,int m) {
    for(int i=0;i<n;i++) {
        vi a;
        for(int j=0;j<m;j++) {
            int t;cin >> t;
            a.push_back(t);
        }
        matrix.push_back(a);
    }
}
/***************************************************************************************************
            function to perform inner product of vectors and compare it with a threshold value
***************************************************************************************************/
vi comp(vii a,vii b,float t_val) {
    vi ret;
    for(int i=0;i<a.size();i++) {
        float sum=0;
        for(int j=0;j<a[i].size();j++) {
            sum+=(a[i][j]*b[i][j]);
        }
        ret.push_back((sum>=t_val?1:0));
    }
    return ret;
}
int main() {
    vii mat,fmat;
    int n,m;
    cin >> n >> m;
    for(int i=0;i<n;i++) {
        vi a=solve(m+1);
        mat.push_back(a);
    }
    cout << endl << "random matrix generated using Blum Blum Shub algorithm:" << endl << endl;
    for(int i=0;i<n;i++) {
        for(int j=0;j<m;j++) printf("%0.3f ",mat[i][j]);
        cout << endl;
    }
    cout << endl;
    mat=orthonormal(mat);
    cout << endl << "Matrix obtained after performing Gram Schmidt Orthonormalisation:" << endl << endl;
    for(int i=0;i<n;i++) {
        //for(int j=0;j<m;j++) cout << mat[i][j] << " ";
        for(int j=0;j<m;j++) printf("%0.4f ",mat[i][j]);
        cout << endl;
    }
    inp(fmat,n,m);
    float threshold;
    cin >> threshold;
    vi ans = comp(mat,fmat,threshold);
    for(int i=0;i<ans.size();i++) cout << ans[i] << ((i==ans.size()-1)?"\n":" ");
    return 0;
}
