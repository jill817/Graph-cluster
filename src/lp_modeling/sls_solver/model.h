#pragma once

#include "solver.h"

/**
 * @brief Data Format Description
 * 1. integer dataset:  user id | supply id | cast times | demand_list (split by ,)
 * 2. demand dataset: demand id | demand amount
 */
namespace solver {
    struct supply_demands {
        int supply_id, demand_id1, demand_id2;
        supply_demands(int _supply, int demand1, int demand2): supply_id(_supply), demand_id1(demand1), demand_id2(demand2)
        {

        }

        bool operator==(supply_demands const & other) const {
            return supply_id == other.supply_id && (
                (demand_id1 == other.demand_id1 && demand_id2 == other.demand_id2) 
                ||
                (demand_id1 == other.demand_id2 && demand_id2 == other.demand_id1)
            );
        }
    };

    using supply_demands_set = std::unordered_set<supply_demands *>;

    class sls_model {
    private:
        opt_solver m_sls_solver;

        std::string integer_path, demand_path;

        std::unordered_map<std::string, int> demand_id_amount; // demand id --> demand amount
        std::unordered_map<std::string, std::unordered_map<std::string, int> > user_supply_cast_times; // user --> supply --> cast times
        std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_set<std::string> > > user_supply_demand_set; // user --> supply --> demand lists
        std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_set<std::string> > > supply_user_demand_set; // supply --> user --> demand lists

        std::unordered_map<std::string, int> var_name_index; // var string --> var index
        std::unordered_map<std::string, std::vector<int> > demand_vars; // demand id --> var index
        std::unordered_map<std::string, std::unordered_set<std::string> > supply_demand_set; // supply --> demand lists
        std::unordered_set<std::string> supply_with_more_demand;

        // id convert
        std::unordered_map<std::string, std::string> supply_id_convert, demand_id_convert, user_id_convert;

        void read_demand_data();
        void read_integer_data();
        void generate_id_convert();
        void model_problem(std::string str, int nia_num, std::string);
    public:
        sls_model(std::string integer_set, std::string demand_set);
        void solve_problem(std::string str, int nia_num, std::string);
        void demo_solve(std::string str);
        void demo_write_lp(std::string str);
    };
};