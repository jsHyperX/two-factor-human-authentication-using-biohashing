/*******************************************************************************************************
                        BASE BIO-HASHING METHOD
*******************************************************************************************************/
#include<bits/stdc++.h>
#define _ ios_base::sync_with_stdio(0);cin.tie(0);
using namespace std;
#define info(x) cerr << #x << " is " << x << endl;
typedef vector<float> vi;
typedef vector<vi> vii;
const int N=100001;
const float pi=3.141592;
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
    for(int i=2;i<(int)sqrt(N);i++) {
        if(prime[i]) {
            pt.push_back(i);            //storing the prime number i in vector pt
            for(int j=i*i;j<N;j+=i) prime[j]=0;
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
vi solve(int m,int seed) {
    vi ans;
    ll n, p[2];
    ll si=pt.size()-1;
    int count=0;
    while(count<2) {
        int a=rand()%si;                // generating a random index in the range of size of array pt
        ll temp=pt[si-a];               // picking up a prime number using the generated index
        if(temp%4==3) { p[count]=temp;count++; }            // checking if the prime obtained is blum or not
    }
    n = ((1LL)*p[0])*p[1];
    for(int i=2;i<=m;i++) {
        ans.push_back(gen(seed,i,p[0],p[1]));
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
                function to perform 1 Dimensional Discrete Fourier transform
**************************************************************************************************/
vi dct(vi seq1) {
    vi seq2;
    int len=seq1.size();
    float factor=sqrt(2.0/(len-1));
    for(int i=0;i<len;i++) {
        float xi=seq1[0];
        xi=((i&1)?xi-seq1[len-1]:xi+seq1[len-1]);
        xi/=2;
        for(int j=1;j<len-1;j++) {
            xi+=seq1[j]*cos((pi*j*i)/len-1);
        }
        xi*=factor;
        seq2.push_back(xi);
    }
    return seq2;
}
/***************************************************************************************************
            function to perform inner product of vectors and compare it with a threshold value
***************************************************************************************************/
vi comp(vii a,vi fmat,float t_val) {
    vi ret;
    for(int i=0;i<a.size();i++) {
        float sum=0;
        for(int j=0;j<a[i].size();j++) {
            sum+=(a[i][j]*fmat[j]);
        }
        //cout << sum << " ";
        ret.push_back((sum>=t_val?1:0));
    }
    cout << endl;
    return ret;
}
/**********************************************************************************************
                    MAIN PROGRAM
**********************************************************************************************/
int main() {
    memset(prime,1,sizeof(prime));              //initializing the array prime with all 1
    gen_p();                                    //generating primes
    vii mat;
    int n,m,seed;
    cout << "Enter the size of the bioemtric feature data:" << endl;
    cin >> n;
    vi tmp,fmat;
    srand (time(NULL));
    for(int i=0;i<n;i++) {
        float x=(float)(rand()%10000);
        tmp.push_back(x);
    }
    cout << endl << "feature matrix before DCT:" << endl << endl;
    for(int i=0;i<n;i++) cout << tmp[i] << " \n"[i==n-1];
    cout << endl;
    fmat=dct(tmp);
    cout << endl << "feature matrix after DCT:" << endl << endl;
    for(int i=0;i<n;i++) cout << fmat[i] << " \n"[i==n-1];
    cout << endl;
    cout << "Enter the number of the random sequence you want to generate:" << endl;
    cin >> m;
    cout << "Enter the Hash Key value:" << endl;
    cin >> seed;
    for(int i=0;i<m;i++) {
        vi a=solve(n+1,seed);
        mat.push_back(a);
    }
    cout << endl << "random matrix generated using Blum Blum Shub algorithm:" << endl << endl;
    for(int i=0;i<m;i++) {
        for(int j=0;j<n;j++) printf("\t%0.4f",mat[i][j]);
        cout << endl;
    }
    cout << endl;
    mat=orthonormal(mat);
    cout << endl << "Matrix obtained after performing Gram Schmidt Orthonormalisation:" << endl << endl;
    for(int i=0;i<m;i++) {
        for(int j=0;j<n;j++) printf("\t%0.4f",mat[i][j]);
        cout << endl;
    }
    cout << endl;
    float threshold;
    cout << "Enter the threshold value :" << endl;
    cin >> threshold;
    vi ans = comp(mat,fmat,threshold);
    cout << "The bio-Hash code of the user is :" << endl;
    for(int i=0;i<ans.size();i++) cout << ans[i] << " \n"[i==ans.size()-1];
    return 0;
}
