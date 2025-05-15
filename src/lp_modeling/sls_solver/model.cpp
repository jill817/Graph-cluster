#include "model.h"
#include "solver.h"

#include <iostream>
#include <fstream>

/**
 * @brief Data Format Description
 * 1. integer dataset:  user id | supply id | cast times | demand_list (split by ,)
 * 2. demand dataset: demand id | demand amount
 * 3. var's name change
 */

namespace solver
{
    void split(std::string &str, std::vector<std::string> &vec, char x)
    {
        vec.clear();
        int start = 0, end = 0;
        for (int i = 0; i < str.length(); i++)
        {
            if (str[i] == x)
            {
                end = i;
                vec.push_back(str.substr(start, end - start));
                start = i + 1;
            }
        }
        end = str.length();
        if (end - start > 0)
        {
            vec.push_back(str.substr(start, end - start));
        }
    }

    sls_model::sls_model(std::string integer_set, std::string demand_set) : integer_path(integer_set), demand_path(demand_set) {}

    void sls_model::read_demand_data()
    {
        demand_id_amount.clear();
        std::cout << "read demand data\n";
        std::ifstream demand_file(demand_path);
        std::string record;
        while (getline(demand_file, record))
        {
            std::vector<std::string> strs;
            split(record, strs, '`');
            demand_id_amount[strs[0]] = std::stoi(strs[1]);
            // demand_id_amount[strs[0]] = std::stoi(strs[1]) * 1.0 / 50000; // 5w original
        }
        std::cout << "demand data size: " << demand_id_amount.size() << std::endl;
        demand_file.close();
    }

    void sls_model::read_integer_data()
    {
        user_supply_cast_times.clear();
        user_supply_demand_set.clear();
        supply_user_demand_set.clear();
        supply_demand_set.clear();
        supply_with_more_demand.clear();
        std::ifstream integer_file(integer_path);
        std::string record;
        // user`supply`demand1,demand2
        while (getline(integer_file, record))
        {
            // only read ?%
            // if (rand() % 10000 < 5000)
            // {
            //     continue;
            // }
            std::vector<std::string> strs;
            split(record, strs, '`');
            user_supply_cast_times[strs[0]][strs[1]] = std::stoi(strs[2]);
            std::vector<std::string> demand_list;
            split(strs[3], demand_list, ',');
            for (auto demand_id : demand_list)
            {
                if (user_supply_demand_set[strs[0]][strs[1]].count(demand_id) == 0)
                {
                    user_supply_demand_set[strs[0]][strs[1]].insert(demand_id);
                }
                if (supply_user_demand_set[strs[1]][strs[0]].count(demand_id) == 0)
                {
                    supply_user_demand_set[strs[1]][strs[0]].insert(demand_id);
                }
                if (supply_demand_set[strs[1]].count(demand_id) == 0)
                {
                    supply_demand_set[strs[1]].insert(demand_id);
                }
            }
        }
        for (auto ele : supply_demand_set)
        {
            if (ele.second.size() >= 2)
            {
                supply_with_more_demand.insert(ele.first);
            }
        }
        integer_file.close();
    }

    /**
     * All ids or strs are used with converted ids
     */
    void sls_model::model_problem(std::string lp_file, int nia_num, std::string result_file)
    {
        std::cout << "model problem\n";
        std::cout << "user num | supply num | demand num\n";
        std::cout << user_supply_cast_times.size() << " | " << supply_demand_set.size() << " | " << demand_id_amount.size() << std::endl;
        
        // rand();//nor for original
        bool enable_rand_query = true;
        // enable_rand_query = false;
        std::string query_str;
        if (enable_rand_query)
        {
            // rand();
            int query_id = rand() % demand_id_convert.size();
            // Random a query demand node
            int cnt = 0;
            for (auto it = demand_id_convert.begin(); it != demand_id_convert.end(); it++)
            {
                if (cnt == query_id)
                {
                    query_str = it->second;
                    break;
                }
                cnt++;
            }
            std::cout << "Using random query selection." << std::endl;
            std::cout << "Selected query string: " << query_str << std::endl;
        }
        else
        {
            // Demand with Maximum var count
            rand();
            std::unordered_map<std::string, int> demand_var_count;
            for (auto ele1 : user_supply_demand_set)
            {
                auto user_id = ele1.first;
                for (auto ele2 : ele1.second)
                {
                    auto supply_id = ele2.first;
                    for (auto demand_id : ele2.second)
                    {
                        demand_var_count[demand_id]++;
                    }
                }
            }
            int max_count = 0;
            std::vector<int> demand_nums;
            for (auto ele : demand_var_count)
            {
                demand_nums.push_back(ele.second);
                if (ele.second > max_count)
                {
                    max_count = ele.second;
                    query_str = ele.first;
                }
            }
            std::cout << demand_nums.size() << std::endl;
            sort(demand_nums.begin(), demand_nums.end());
            reverse(demand_nums.begin(), demand_nums.end());
            int top_n_count = demand_nums[10];
            std::cout << top_n_count << std::endl;
            for (auto ele : demand_var_count)
            {
                if (ele.second == top_n_count)
                {
                    query_str = ele.first;
                }
            }
            std::cout << "Using variable count-based query selection." << std::endl;
            std::cout << "Selected query string: " << query_str << std::endl;
        }

        var_name_index.clear();
        demand_vars.clear();
        // register vars
        for (auto ele1 : user_supply_demand_set)
        {
            auto user_id = ele1.first;
            for (auto ele2 : ele1.second)
            {
                auto supply_id = ele2.first;
                for (auto demand_id : ele2.second)
                {
                    std::string var_name = "x_" + user_id_convert[user_id] + "_" + supply_id_convert[supply_id] + "_" + demand_id_convert[demand_id];
                    double var_coeff;
                    if (demand_id_convert[demand_id] == query_str)
                    {
                        var_coeff = -1.0;
                    }
                    else
                    {
                        var_coeff = 0;
                    }
                    int var_index = m_sls_solver.register_var(var_name, var_coeff, true, true, 0, user_supply_cast_times[user_id][supply_id]);
                    var_name_index[var_name] = var_index;
                    demand_vars[demand_id].push_back(var_index);
                }
            }
        }
        std::unordered_map<std::string, int> slack_vars;
        double BIG_COEFFICIENT = 31708040;
        for (auto ele : demand_vars) {
            auto demand_id = ele.first;
            std::string slack_var_name = "s_" + demand_id_convert[demand_id];
            int slack_var_index = m_sls_solver.register_var(slack_var_name, BIG_COEFFICIENT, true, true, 0, INFINITY);
            slack_vars[demand_id] = slack_var_index;
        }
        

        std::cout << "register vars done\n";

        // Part I. Basic Constraints: Linear
        // 2.1 user-supply: sum of vars <= cast count
        for (auto ele1 : user_supply_demand_set)
        {
            auto user_id = ele1.first;
            for (auto ele2 : ele1.second)
            {
                auto supply_id = ele2.first;
                int cast_count = user_supply_cast_times[user_id][supply_id];
                // if(cast_count>50) {
                //     std::cout<<user_id<<std::endl;
                //     std::cout<<supply_id<<std::endl;
                //     std::cout<<" >50 "<<std::endl;
                // }
                monomial_vector m_monomials;
                for (auto demand_id : ele2.second)
                {
                    std::string var_name = "x_" + user_id_convert[user_id] + "_" + supply_id_convert[supply_id] + "_" + demand_id_convert[demand_id];
                    int_table vars;
                    vars.insert(var_name_index[var_name]);
                    m_monomials.push_back(monomial(1.0, vars));
                }
                polynomial poly(m_monomials);
                m_sls_solver.register_constraint(poly, constraint_kind::CAST, "user_supply_" + user_id + "_" + supply_id, false, true, 0, cast_count);
            }
        }
        std::cout << "2.1 supply done\n";

        // // 2.2 demand: sum of vars >= demand need
        // for (auto ele : demand_vars)
        // {
        //     monomial_vector m_monomials;
        //     auto demand_id = ele.first;
        //     if (demand_id == query_str)
        //     {
        //         continue;
        //     }
        //     int demand_amount = demand_id_amount[demand_id];
        //     for (auto m_var : ele.second)
        //     {
        //         int_table vars;
        //         vars.insert(m_var);
        //         m_monomials.push_back(monomial(1.0, vars));
        //     }
        //     polynomial poly(m_monomials);
        //     m_sls_solver.register_constraint(poly, constraint_kind::DEMAND, "demand_" + demand_id, true, false, demand_amount, 0);
        // }

        // 2.2 demand: sum of vars + slack var >= demand need
        for (auto ele : demand_vars) {
            monomial_vector m_monomials;
            auto demand_id = ele.first;
            if (demand_id == query_str) {
                std::cout << "Skipping demand " << demand_id << std::endl;
                continue;
            }
            int demand_amount = demand_id_amount[demand_id];
        
            // 添加原有变量
            for (auto m_var : ele.second) {
                int_table vars;
                vars.insert(m_var);
                m_monomials.push_back(monomial(1.0, vars));
            }
        
            // 添加辅助变量
            int_table slack_vars_table;
            slack_vars_table.insert(slack_vars[demand_id]);
            m_monomials.push_back(monomial(1.0, slack_vars_table));
        
            polynomial poly(m_monomials);
            m_sls_solver.register_constraint(poly, constraint_kind::DEMAND, "demand_" + demand_id, true, false, demand_amount, 0);
        }
        
        std::cout << "2.2 demand done\n";

        // Part II. Advanced Constraints: Multilinear
        // 3.7 random choose a supply node and two demand nodes
        // int num_3_7 = nia_num;
        // supply_demands_set m_supply_demands;
        // for (int k = 0; k < num_3_7; k++)
        // {
        //     std::string supply_str, demand_str1, demand_str2;
        //     int supply_id, demand_id1, demand_id2;
        //     while (true)
        //     { // in case we do not random to a right node
        //         // 1. random choose a supply node
        //         int supply_size = supply_with_more_demand.size();
        //         supply_id = rand() % supply_size;
        //         int cnt = 0;
        //         for (auto it = supply_with_more_demand.begin(); it != supply_with_more_demand.end(); it++)
        //         {
        //             if (cnt == supply_id)
        //             {
        //                 supply_str = *it;
        //             }
        //             cnt++;
        //         }
        //         ASSERT(supply_demand_set[supply_str].size() >= 2);
        //         // 2. random choose two demand nodes attached to this supply node
        //         do
        //         {
        //             demand_id1 = rand() % supply_demand_set[supply_str].size();
        //             demand_id2 = rand() % supply_demand_set[supply_str].size();
        //         } while (demand_id1 == demand_id2);
        //         cnt = 0;
        //         supply_demands *curr_supply_demand = new supply_demands(supply_id, demand_id1, demand_id2);
        //         if (m_supply_demands.count(curr_supply_demand) != 0)
        //         {
        //             continue;
        //         }
        //         m_supply_demands.insert(curr_supply_demand);

        //         for (auto it = supply_demand_set[supply_str].begin(); it != supply_demand_set[supply_str].end(); it++)
        //         {
        //             if (cnt == demand_id1)
        //             {
        //                 demand_str1 = *it;
        //             }
        //             if (cnt == demand_id2)
        //             {
        //                 demand_str2 = *it;
        //             }
        //             cnt++;
        //         }
        //         break;
        //     }
        //     std::cout << "3.7 supply id | 3.7 demand id1 | 3.7 demand id2\n";
        //     std::cout << supply_str << " | " << demand_str1 << " | " << demand_str2 << std::endl;
        //     // 3. generate constraint
        //     /*
        //         sum_i x{ui, s, d1} / sum_ij x{ui, sj, d1} >= sum_i x{ui, s, d2} / sum_ij x{ui, sj, d2}
        //         -------------------------------------------------------------------------------------
        //         sum_i x{ui, s, d1} * sum_ij x{ui, sj, d2} >= sum_i x{ui, s, d2} * sum_ij x{ui, sj, d1}
        //         ------------------   --------------------    ------------------   --------------------
        //             p1                     p4                     p2                    p3
        //     */
        //     // loop all edges
        //     monomial_vector p1s, p2s, p3s, p4s;
        //     int_table m_vars;
        //     for (auto ele1 : user_supply_demand_set)
        //     {
        //         auto curr_user_id = ele1.first;
        //         for (auto ele2 : ele1.second)
        //         {
        //             auto curr_supply_id = ele2.first;
        //             for (auto curr_demand_id : ele2.second)
        //             {
        //                 if (curr_supply_id == supply_str && curr_demand_id == demand_str1)
        //                 { // p1 case
        //                     m_vars.clear();
        //                     m_vars.insert(var_name_index["x_" + user_id_convert[curr_user_id] + "_" + supply_id_convert[curr_supply_id] + "_" + demand_id_convert[curr_demand_id]]);
        //                     p1s.push_back(monomial(1.0, m_vars));
        //                 }
        //                 if (curr_supply_id == supply_str && curr_demand_id == demand_str2)
        //                 { // p2 case
        //                     m_vars.clear();
        //                     m_vars.insert(var_name_index["x_" + user_id_convert[curr_user_id] + "_" + supply_id_convert[curr_supply_id] + "_" + demand_id_convert[curr_demand_id]]);
        //                     p2s.push_back(monomial(1.0, m_vars));
        //                 }
        //                 if (curr_demand_id == demand_str1)
        //                 { // p3 case
        //                     m_vars.clear();
        //                     m_vars.insert(var_name_index["x_" + user_id_convert[curr_user_id] + "_" + supply_id_convert[curr_supply_id] + "_" + demand_id_convert[curr_demand_id]]);
        //                     p3s.push_back(monomial(1.0, m_vars));
        //                 }
        //                 if (curr_demand_id == demand_str2)
        //                 { // p4 case
        //                     m_vars.clear();
        //                     m_vars.insert(var_name_index["x_" + user_id_convert[curr_user_id] + "_" + supply_id_convert[curr_supply_id] + "_" + demand_id_convert[curr_demand_id]]);
        //                     p4s.push_back(monomial(1.0, m_vars));
        //                 }
        //             }
        //         }
        //     }
        //     polynomial p1(p1s), p2(p2s), p3(p3s), p4(p4s);
        //     // std::cout << "four polys' size: " << p1.size() << ", " << p2.size() << ", " << p3.size() << ", " << p4.size() << std::endl;
        //     p1.mul(p4);
        //     p2.mul(p3);
        //     // std::cout << "left and right poly's size: " << p1.size() << ", " << p2.size() << std::endl;
        //     polynomial cons_poly = p1 - p2;
        //     // std::cout << "cons size: " << cons_poly.size() << std::endl;
        //     m_sls_solver.register_constraint(cons_poly, constraint_kind::DEMAND_COMPARE, "cons_3_7" + supply_str + "_" + demand_str1 + "_" + demand_str2, true, 0, false, 0);
        // }

        std::cout << "start solve\n";
        m_sls_solver.write_lp_file(lp_file);
        std::cout << "write lp file done\n";
        exit(0);
        // m_sls_solver.set_cutoff(300);
        // m_sls_solver.solve();
        // m_sls_solver.display_results(std::cout);
    }

    void sls_model::generate_id_convert()
    {
        user_id_convert.clear();
        supply_id_convert.clear();
        demand_id_convert.clear();
        int idx;
        for (auto ele1 : user_supply_demand_set)
        {
            auto user_id = ele1.first;
            if (user_id_convert.count(user_id) == 0)
            {
                idx = user_id_convert.size();
                user_id_convert[user_id] = std::to_string(idx);
            }
            for (auto ele2 : ele1.second)
            {
                auto supply_id = ele2.first;
                if (supply_id_convert.count(supply_id) == 0)
                {
                    idx = supply_id_convert.size();
                    supply_id_convert[supply_id] = std::to_string(idx);
                }
                for (auto demand_id : ele2.second)
                {
                    if (demand_id_convert.count(demand_id) == 0)
                    {
                        idx = demand_id_convert.size();
                        demand_id_convert[demand_id] = std::to_string(idx);
                    }
                }
            }
        }
    }

    void sls_model::solve_problem(std::string lp_file, int nia_num, std::string result_file)
    {
        std::cout << "start read demand data\n";
        read_demand_data();
        std::cout << "read demand data done\n";
        std::cout << "start read supply data\n";
        read_integer_data();
        std::cout << "read supply data done\n";
        std::cout << "start generate id convert\n";
        generate_id_convert();
        std::cout << "generate id convert done\n";
        std::cout << "start model problem\n";
        model_problem(lp_file, nia_num, result_file);
        std::cout << "model problem done\n";
    }

    /**
     * @brief Demo for write lp file
     *
     * Minimize -3 x1 - 5 x2 - 2 x3 + x4
     *
     * x1 + 2 x2 * x3 + x2 <= 15
     * 2 x1 + 4 x2 + 5 x3 + x4 >= 10
     * x1 - 3 x2 + 4 x3 + x4 = 5
     * x2 + [3 x2 * x3  - 4 x1 * x2 ] <= 10
     *
     * x1 in [0, 10]
     * x2 in [0, 5]
     * x3 <= 3
     * x4 >= 0s
     */
    void sls_model::demo_solve(std::string str)
    {
        m_sls_solver.register_var("x1", -3, true, true, 0, 10);
        m_sls_solver.register_var("x2", -5, true, true, 0, 5);
        m_sls_solver.register_var("x3", -2, false, true, 0, 3);
        m_sls_solver.register_var("x4", 1, true, false, 0, 0);

        long_double_vector coeff_vec = {1.0, 2.0, 1.0};
        int_table_vector mono_vars = {{0}, {1, 2}, {1}};
        polynomial poly1(coeff_vec, mono_vars);
        m_sls_solver.register_constraint(poly1, constraint_kind::OTHER, "c1", false, true, 0, 15);
        coeff_vec = {2.0, 4.0, 5.0, 1.0};
        mono_vars = {{0}, {1}, {2}, {3}};
        polynomial poly2(coeff_vec, mono_vars);
        m_sls_solver.register_constraint(poly2, constraint_kind::OTHER, "c2", true, false, 10, 0);
        coeff_vec = {1.0, -3.0, 4.0, 1.0};
        mono_vars = {{0}, {1}, {2}, {3}};
        polynomial poly3(coeff_vec, mono_vars);
        m_sls_solver.register_constraint(poly3, constraint_kind::OTHER, "c3", true, true, 5, 5);
        coeff_vec = {1.0, 3.0, -4.0};
        mono_vars = {{1}, {1, 2}, {0, 1}};
        polynomial poly4(coeff_vec, mono_vars);
        m_sls_solver.register_constraint(poly4, constraint_kind::OTHER, "c4", false, true, 0, 10);

        m_sls_solver.set_cutoff(0.5);
        // m_sls_solver.set_stepoff(1000);
        m_sls_solver.solve();
        m_sls_solver.display_results(std::cout);
    }

    /**
     * @brief Demo for write lp file
     *
     * Minimize -3 x1 - 5 x2 - 2 x3 + x4
     *
     * x1 + 2 x2 * x3 + x2 <= 15
     * 2 x1 + 4 x2 + 5 x3 + x4 >= 10
     * x1 - 3 x2 + 4 x3 + x4 = 5
     * x2 + [3 x2 * x3  - 4 x1 * x2 ] <= 10
     *
     * x1 in [0, 10]
     * x2 in [0, 5]
     * x3 <= 3
     * x4 >= 0
     */
    void sls_model::demo_write_lp(std::string str)
    {
        m_sls_solver.register_var("x1", -3, true, true, 0, 10);
        m_sls_solver.register_var("x2", -5, true, true, 0, 5);
        m_sls_solver.register_var("x3", -2, false, true, 0, 3);
        m_sls_solver.register_var("x4", 1, true, false, 0, 0);

        long_double_vector coeff_vec = {1.0, 2.0, 1.0};
        int_table_vector mono_vars = {{0}, {1, 2}, {1}};
        polynomial poly1(coeff_vec, mono_vars);
        m_sls_solver.register_constraint(poly1, constraint_kind::OTHER, "c1", false, true, 0, 15);
        coeff_vec = {2.0, 4.0, 5.0, 1.0};
        mono_vars = {{0}, {1}, {2}, {3}};
        polynomial poly2(coeff_vec, mono_vars);
        m_sls_solver.register_constraint(poly2, constraint_kind::OTHER, "c2", true, false, 10, 0);
        coeff_vec = {1.0, -3.0, 4.0, 1.0};
        mono_vars = {{0}, {1}, {2}, {3}};
        polynomial poly3(coeff_vec, mono_vars);
        m_sls_solver.register_constraint(poly3, constraint_kind::OTHER, "c3", true, true, 5, 5);
        coeff_vec = {1.0, 3.0, -4.0};
        mono_vars = {{1}, {1, 2}, {0, 1}};
        polynomial poly4(coeff_vec, mono_vars);
        m_sls_solver.register_constraint(poly4, constraint_kind::OTHER, "c4", false, true, 0, 10);

        m_sls_solver.write_lp_file(str);
    }
};