#include<bits/stdc++.h>
#define _ ios_base::sync_with_stdio(0);cin.tie(0);
using namespace std;
typedef vector<double> vi;
typedef vector<vi> vii;
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
vii orthogonal(vii vec) {
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
int main() {
    vii mat;
    int n,m;
    cin >> n >> m;
    for(int i=0;i<n;i++) {
        vi a;
        for(int j=0;j<m;j++) {
            int t;cin >> t;
            a.push_back(t);
        }
        mat.push_back(a);
    }
    mat=orthogonal(mat);
    for(int i=0;i<mat.size();i++) {
        for(int j=0;j<mat[i].size();j++) cout << mat[i][j] << " ";
        cout << endl;
    }
    return 0;
}
