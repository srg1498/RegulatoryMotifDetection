#include<bits/stdc++.h>
#define N 100

using namespace std;

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

pair<int,int> distance(vector<int> s,string dna,int n,int l,int flag)
{
    int hd = INT_MAX;
    int st;
    pair<int,int> p;

    if(flag == 1)
    {
        for(int i=0;i<n-l+1;i++)
        {
            int mismatches = 0;
            for(int j=1;j<=l;j++)
            {
                if(s[j] == 1)
                {
                    if((dna[i+j-1]!='A'))
                        mismatches++;
                }
                else if(s[j] == 2)
                {
                    if((dna[i+j-1]!='C'))
                        mismatches++;
                }
                else if(s[j] == 3)
                {
                    if((dna[i+j-1]!='G'))
                        mismatches++;
                }
                else
                {
                    if((dna[i+j-1]!='T'))
                        mismatches++;
                }
            }
            if(mismatches < hd)
            {
                hd = mismatches;
                st = i+1;
            }
        }
    }
    else
    {
        for(int i=0;i<n-l+1;i++)
        {
            int mismatches = 0;
            for(int j=1;j<=l;j++)
            {
                if(s[j] == 1)
                {
                    if((dna[i+j-1]!='a'))
                        mismatches++;
                }
                else if(s[j] == 2)
                {
                    if((dna[i+j-1]!='c'))
                        mismatches++;
                }
                else if(s[j] == 3)
                {
                    if((dna[i+j-1]!='g'))
                        mismatches++;
                }
                else
                {
                    if((dna[i+j-1]!='t'))
                        mismatches++;
                }
            }
            if(mismatches < hd)
            {
                hd = mismatches;
                st = i+1;
            }
        }
    }
    p = make_pair(hd,st);

    return p;
}

pair< int,vector<int> > total_distance(vector<int> s,string DNA[N],int t,int n,int l,int flag)
{
    int dist = 0;
    pair< int,vector<int> > p;
    pair<int,int> q;

    for(int i=1;i<=t;i++)
    {
        q = distance(s,DNA[i],n,l,flag);
        dist = dist + q.first;
        p.second.push_back(q.second);
    }
    p.first = dist;

    return p;
}

pair< string,vector<int> > median_string_motif_search(string DNA[N],int t,int n,int l,int flag)
{
    int best_distance;
    string best_motif;
    vector<int> best_motif_int(l+1);
    vector<int> s(l+1);
    vector<int> z(l+1);
    vector<int> all_motifs(t);
    pair< int,vector<int> > p;
    pair< string,vector<int> > final_result;

    for(int i=0;i<=t;i++)
        s[i] = z[i] = 1;
    p = total_distance(s,DNA,t,n,l,flag);
    best_distance = p.first;
    all_motifs = p.second;
    best_motif_int = s;
    while(true)
    {
        s = next_leaf(s,l,4);
        pair< int,vector<int> > q;
        q = total_distance(s,DNA,t,n,l,flag);
        if(q.first < best_distance)
        {
            best_distance = q.first;
            all_motifs = q.second;
            best_motif_int = s;
        }
        if(s == z)
            break;
    }
    for(int i=1;i<=l;i++)
    {
        if(best_motif_int[i] == 1)
            best_motif.push_back('a');
        else if(best_motif_int[i] == 2)
            best_motif.push_back('c');
        else if(best_motif_int[i] == 3)
            best_motif.push_back('g');
        else
            best_motif.push_back('t');
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
    for(int i=0;i<t;i++)
        cout<<result.second[i]<<" ";
}

int main()
{
    int t,n,l;
    int flag;
    char ques;
    string DNA[N];
    pair< string,vector<int> > result;

    cout<<"Enter total number of dna sequences: ";
    cin>>t;
    cout<<"\nEnter the length of the dna sequence: ";
    cin>>n;
    cout<<"\nEnter the length of the motif: ";
    cin>>l;
    cout<<"\nDo you want the dna sequences in capital letters(y/n): ";
    cin>>ques;
    flag = (ques=='y') ? 1 : 0;
    cout<<"\nEnter all the dna sequences:\n\n";
    read_all_dna_sequences(DNA,t);

    result = median_string_motif_search(DNA,t,n,l,flag);
    final_result(result,t);
    cout<<"\n\n\n";

    return 0;
}
