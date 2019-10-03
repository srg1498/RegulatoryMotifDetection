#include<bits/stdc++.h>
#define N 100

using namespace std;

char nucleotides_char[] = {'a','c','g','t'};

vector<int> next_leaf(vector<int> s,int L,int k)
{
    for(int i=L;i>=1;i--)
    {
        if(s[i] < k)
        {
            s[i] = s[i] + 1;
            return s;
        }
        s[i] = 1;
    }
    return s;
}

int count_of_character(char c,int j,vector<int> s,string DNA[N],int t)
{
    int cnt = 0;
    for(int i=1;i<=t;i++)
    {
        if((DNA[i][s[i]+j-2]==c) || ((DNA[i][s[i]+j-2]==c-32)))
            cnt++;
    }
    return cnt;
}

pair<string,int> score(vector<int> s,string DNA[N],int t,int l)
{
    int sc=0;
    string consensus_motif;
    pair<string,int> p;

    for(int col=1;col<=l;col++)
    {
        int max_count = INT_MIN;
        char max_char;
        for(int c=0;c<4;c++)
        {
            int cnt = count_of_character(nucleotides_char[c],col,s,DNA,t);
            if(cnt > max_count)
            {
                max_count = cnt;
                max_char = nucleotides_char[c];
            }
        }
        sc += max_count;
        consensus_motif.push_back(max_char);
    }
    p = make_pair(consensus_motif,sc);

    return p;
}

pair< string,vector<int> > brute_force_motif_search(string DNA[N],int t,int n,int l)
{
    int best_score;
    string best_motif;
    vector<int> s(t+1);
    vector<int> z(t+1);
    vector<int> all_motifs(t+1);
    pair<string,int> p;
    pair< string,vector<int> > final_result;

    for(int i=0;i<=t;i++)
        s[i] = z[i] = 1;
    p = score(s,DNA,t,l);
    best_motif = p.first;
    best_score = p.second;
    all_motifs = s;
    while(true)
    {
        s = next_leaf(s,t,n-l+1);
        pair<string,int> q;
        q = score(s,DNA,t,l);
        if(q.second > best_score)
        {
            best_score = q.second;
            best_motif = q.first;
            all_motifs = s;
        }
        if(s == z)
            break;
    }
    final_result = make_pair(best_motif,all_motifs);

    return final_result;
}

void read_all_dna_sequences(string DNA[N],int t)
{
    for(int i=1;i<=t;i++)
    {
        cout<<"Sequence number "<<i<<": ";
        cin>>DNA[i];
        cout<<"\n";
    }
    cout<<"\n\n";
}

void final_result(pair< string,vector<int> > result,int t)
{
    cout<<"Final motif: "<<result.first<<"\n";
    cout<<"All the starting positions of each motif in the dna sequences: ";
    for(int i=1;i<=t;i++)
        cout<<result.second[i]<<" ";
}

int main()
{
    int t,n,l;
    string DNA[N];
    pair< string,vector<int> > result;

    cout<<"Enter total number of dna sequences: ";
    cin>>t;
    cout<<"\nEnter the length of the dna sequence: ";
    cin>>n;
    cout<<"\nEnter the length of the motif: ";
    cin>>l;
    cout<<"\nEnter all the dna sequences:\n\n";
    read_all_dna_sequences(DNA,t);

    result = brute_force_motif_search(DNA,t,n,l);
    final_result(result,t);
    cout<<"\n\n\n";

    return 0;
}
