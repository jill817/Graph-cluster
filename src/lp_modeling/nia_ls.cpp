#include "sls_solver/model.h"

std::string date[5] = {"405", "505", "515", "527", "000"};
const int date_index = 4;
// std::string demand_path = "data_by_date/demand/demand_20230" + date[date_index];
// std::string integer_path = "data_by_date/sample/integer_sample_100000_20230" + date[date_index];

std::string demand_path = "/pub/netdisk1/lijy/alimama/data/demand_330_0416.txt";
std::string integer_path = "/pub/netdisk1/lijy/alimama/data/supply_0416_cluster6_format.txt";

int nia_num = 100;

/*
    Minimize -3 x1 - 5 x2 - 2 x3 + x4

    x1 + 2 x2 * x3 + x2 <= 15
    2 x1 + 4 x2 + 5 x3 + x4 >= 10
    x1 - 3 x2 + 4 x3 + x4 = 5

    x1 in [0, 10]
    x2 in [0, 5]
    x3 <= 3
    x4 >= 0
*/

void test_write_lp()
{
    solver::sls_model m_model(integer_path, demand_path);
    std::string str = "demo_lp.lp";
    m_model.demo_write_lp(str);
}

void test_solve_demo()
{
    solver::sls_model m_model(integer_path, demand_path);
    m_model.demo_solve("demo_solve.lp");
}

void solve()
{
    solver::sls_model m_model(integer_path, demand_path);
    std::cout<<"start solve"<<std::endl;
    m_model.solve_problem("/pub/netdisk1/lijy/alimama/lp_learning/0416_connect_cluster6.lp", 
                            nia_num, 
                            "experiment/" + std::to_string(nia_num) + ".res");
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        UNREACHABLE();
        return -1;
    }
    if (std::strcmp(argv[1], "solve") == 0)
    {
        solve();
    }
    else if (std::strcmp(argv[1], "dsolve") == 0)
    {
        test_solve_demo();
    }
    else if (std::strcmp(argv[1], "dlp") == 0)
    {
        test_write_lp();
    }
    else if (std::strcmp(argv[1], "exp") == 0)
    {
        nia_num = std::stoi(argv[2]);
        solve();
    }
    return 0;
}
