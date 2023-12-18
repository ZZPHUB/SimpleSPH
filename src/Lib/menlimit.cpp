#include "Lib.H"
using namespace std;

void get_max_men(void)
{
    struct rlimit rs;
    rs.rlim_cur = 0;
    rs.rlim_max = 0;
    getrlimit(RLIMIT_AS,&rs);
    cout << rs.rlim_cur/8*pow(2,30) << "MB " << rs.rlim_max/8*pow(2,30) << "MB" << endl;
}

 
void set_max_men(void)
{
    struct rlimit rs;
    rs.rlim_cur = (rlim_t)pow(2,30);
    rs.rlim_max = (rlim_t)pow(2,30);
 
    setrlimit(RLIMIT_AS, &rs);
    return;
}