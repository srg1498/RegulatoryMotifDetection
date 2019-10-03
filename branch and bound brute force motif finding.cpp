#include<bits/stdc++.h>
#define N 100

using namespace std;

char nucleotides_char[] = {'a','c','g','t'};

pair< vector<int>,int > next_vertex(vector<int> s,int i,int L,int k)
{
    pair< vector<int>,int > p;

    if(i < L)
    {
        s[i+1] = 1;
        p = make_pair(s,i+1);
        return p;
    }
    else
    {
        for(int j=L;j>=1;j--)
        {
            if(s[j] < k)
            {
                s[j] = s[j] + 1;
                p = make_pair(s,j);
                return p;
            }
        }
    }
    p = make_pair(s,0);

    return p;
}

pair< vector<int>,int > bypass(vector<int> s,int i,int k)
{
    pair< vector<int>,int > p;

    for(int j=i;j>=1;j--)
    {
        if(s[j] < k)
        {
            s[j] = s[j] + 1;
            p = make_pair(s,j);
        }
    }
    p = make_pair(s,0);

    return p;
}

int count_of_character(char c,int j,vector<int> s,string DNA[N],int co)
{
    int cnt = 0;
    for(int i=1;i<=co;i++)
    {
        if((DNA[i][s[i]+j-2]==c) || ((DNA[i][s[i]+j-2]==c-32)))
            cnt++;
    }
    return cnt;
}

pair<string,int> score(vector<int> s,string DNA[N],int i,int l)
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
            int cnt = count_of_character(nucleotides_char[c],col,s,DNA,i);
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

pair< string,vector<int> > branch_and_bound_brute_force_motif_search(string DNA[N],int t,int n,int l)
{
    int best_score = 0;
    string best_motif;
    vector<int> all_motifs(t+1);
    vector<int> s(t+1);
    pair<string,int> p;
    pair< string,vector<int> > final_result;

    for(int i=0;i<=t;i++)
        s[i] = 1;
    int i = 1;
    while(i>0)
    {
        if(i < t)
        {
            p = score(s,DNA,i,l);
            int optimistic_score = p.second + (t-i)*l;
            if(optimistic_score < best_score)
            {
                pair< vector<int>,int > q;
                q = bypass(s,i,n-l+1);
                s = q.first;
                i = q.second;
            }
            else
            {
                pair< vector<int>,int > q;
                q = next_vertex(s,i,t,n-l+1);
                s = q.first;
                i = q.second;
            }
        }
        else
        {
            p = score(s,DNA,t,l);
            if(p.second > best_score)
            {
                best_score = p.second;
                best_motif = p.first;
                all_motifs = s;
            }
            pair< vector<int>,int > q;
            q = next_vertex(s,i,t,n-l+1);
            s = q.first;
            i = q.second;
        }
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

    result = branch_and_bound_brute_force_motif_search(DNA,t,n,l);
    final_result(result,t);
    cout<<"\n\n\n";

    return 0;
}

