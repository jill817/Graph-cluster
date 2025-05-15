#include "solver.h"

#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <exception>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <chrono>
#include <sys/types.h>
#include <string>

/**
 * @todo
 * 1. initialize greedy: count first, demand amount next
 * 2. adjust factor polynomial's value when inserting new var
 * 3. insert operation, floor or ceil
 * 4. in case > <, insert value
 *
 * @todo
 * 1. object weight = 0 before find sat, set 1 when find sat
 * 2. store var-constraint
 *
 * @date 23/05/09
 * @todo (code review)
 * 1. check solution sat
 * 2. lp file write object function, space before -1
 * 3. unbounded speed: BMS or other algorithms
 * 4. read lp format function
 *
 * @date 23/05/26
 * @todo (code review) VND
 * 1. for each unsat constraint, insert n-dimensional operation
 * 2. two level: first select var-pair which do not factor each other, second select that factor each other
 *
 * @date 23/05/30
 * @todo (code review) VND Lift
 * 1. random select 10-vnd paired lift move, paired lift move selects vars from objective and do not factor each other in any constraint
 *
 */

namespace solver
{
    struct opt_solver::imp
    {
        // basic data structures
        var_vector m_vars;
        std::unordered_map<std::string, int> m_var_map;
        constraint_vector m_constraints;
        long double m_object_sum;
        long double m_object_weight;
        long double m_sum_weight;
        long double m_avg_weight;
        const long double weight_gap = 1e3;
        assignment m_assignment;
        operation_table m_operation_table;
        int_table m_var_bound_cons, m_cast_cons, m_demand_cons, m_demand_compare_cons, m_other_cons;
        std::unordered_map<int, long double> cast_cons_left;

        // statistics
        int m_lift_step;
        int m_lift_pair_step;
        int m_vnd_lift_step;
        bool find_sat_status = false;
        long double m_best_object_sum = 1e100;
        assignment m_best_obj_assignment;
        int m_best_unsat_num = INT32_MAX;
        assignment m_best_unsat_assignment;
        const int random_walk_num = 3;
        int no_improve_unsat; // no improve count if not find sat
        int no_improve_sat;   // no improve count if find sat
        int m_restart;
        int_table m_vars_in_obj;

        // status
        bool is_sat_status;
        int_table m_unsat_constraints;
        int_table m_unbounded_constraints;

        // control
        int m_step;
        const int max_step = INT32_MAX;
        const int BMS_VAL = 1000;
        const int VND_BMS_VAL = 100;
        const int DOUBLE_BMS_VAL = 50;

        // switch
        const bool enable_bound_constraints = true;
        bool enable_initialize_greedy = false;
        const bool enable_tabu = true;
        const int tabu_const = 3, tabu_rand = 10;
        const bool enable_lift_move = true;
        const bool enable_paired_lift_move = false;
        const bool enable_vnd_lift_move = false;
        const bool enable_all_lift_move_pool = false;

        // VND
        const bool enable_vnd_search = false;
        const bool enable_vnd_two_level = false;
        // gradient_vector                                 m_gradient_vector;
        int gradient_step;
        int m_vnd_step;
        int m_vnd_step_two_level;
        double learning_rate = 1;
        const double laerning_rate_decay = 0.95;

        // time
        bool enable_cutoff = false;
        long double m_cutoff;
        std::chrono::steady_clock::time_point m_start_time;
        bool enable_stepoff = false;
        int m_stepoff;

        // trace
        std::ofstream tout;

        imp() : tout("/Users/zhonghanwang/Desktop/multilinear_solver/ls.trace")
        {
            TRACE(tout << "begin\n";);
        }

        // functions
        /**
         * @brief Register a new var in sls solver
         *
         * @param str name of var
         * @param _coeff coeffcient of this var in objective
         * @param lower has lower bound or not
         * @param upper has upper bound or not
         * @param lb lower bound
         * @param ub upper bound
         * @return int var's index in sls solver
         */
        int register_var(std::string str, long double _coeff, bool lower, bool upper, long double lb, long double ub)
        {
            int index = m_vars.size();
            m_vars.push_back(var(index, str, _coeff, lower, upper, lb, ub));
            m_var_map[str] = index;
            if (_coeff != 0)
            {
                m_vars_in_obj.insert(index);
            }
            return index;
        }

        /**
         * @brief Register a new constraint in sls solver
         */
        void register_constraint(polynomial const &poly, constraint_kind _kind, std::string m_name, bool lower, bool upper, long double lb, long double ub)
        {
            int index = m_constraints.size();
            m_constraints.push_back(new constraint(index, m_name, _kind, poly, lower, upper, lb, ub));
            switch (_kind)
            {
            case constraint_kind::CAST:
                m_cast_cons.insert(index);
                break;
            case constraint_kind::DEMAND:
                m_demand_cons.insert(index);
                break;
            case constraint_kind::DEMAND_COMPARE:
                m_demand_compare_cons.insert(index);
                break;
            case constraint_kind::BOUND:
                m_var_bound_cons.insert(index);
                break;
            case constraint_kind::OTHER:
                m_other_cons.insert(index);
                break;
            default:
                break;
            }
        }

        void set_cutoff(long double b)
        {
            enable_cutoff = true;
            m_cutoff = b;
        }

        void set_stepoff(int _step)
        {
            enable_stepoff = true;
            m_stepoff = _step;
        }

        void enable_greedy_initialize()
        {
            enable_initialize_greedy = true;
        }

        long double time_elapsed() const
        {
            auto m_curr_time = std::chrono::steady_clock::now();
            std::chrono::duration<long double> duration = m_curr_time - m_start_time;
            return duration.count();
        }

        void enable_switch_choices()
        {
            if (enable_bound_constraints)
            {
                add_var_bound_constraints();
            }
        }

        /**
         * @brief Add var's bound as a new constraint
         *
         */
        void add_var_bound_constraints()
        {
            for (auto v : m_vars)
            {
                if (!v.m_lower && !v.m_upper)
                {
                    continue;
                }
                monomial_vector m_monomials;
                int_table vars = {v.m_index};
                m_monomials.push_back(monomial(1.0, vars));
                polynomial poly(m_monomials);
                register_constraint(poly, constraint_kind::BOUND, v.m_name + "_bound", v.m_lower, v.m_upper, v.m_lb, v.m_ub);
            }
        }

        void initialize_no_improve()
        {
            no_improve_sat = 0;
            no_improve_unsat = 0;
        }

        void initialize()
        {
            std::cout << "begin initialize\n";
            m_start_time = std::chrono::steady_clock::now();
            initialize_no_improve();
            initialize_assignment();
            initialize_constraints();
            initialize_var_constraints();
            initialize_objective();
            update_best_information();
            gradient_step = 0;
            m_vnd_step = 0;
            m_vnd_step_two_level = 0;
            std::cout << "end initialize\n";
        }

        /**
         * @brief Check whether sat constraints this step are actually sat
         *  used for debug
         */
        void check_sat_constraints() const
        {
            std::cout << "enter check sat constraints\n";
            for (int i = 0; i < m_constraints.size(); i++)
            {
                if (m_unsat_constraints.count(i) != 0)
                {
                    continue;
                }
                long double sum = m_constraints[i]->m_poly.get_sum();
                auto res = check_constraint_sat_using_sum(m_constraints[i], sum);
                if (!(res.first && res.second))
                {
                    UNREACHABLE();
                    std::cout << m_constraints[i]->m_name << std::endl;
                    std::cout << sum << ", " << m_constraints[i]->m_lower_bound << ", " << m_constraints[i]->m_upper_bound << std::endl;
                    break;
                }
            }
        }

        assignment m_temp_assignment;

        void check_result()
        {
            m_temp_assignment.copy(m_assignment);
            if (find_sat_status)
            {
                check_result_sat();
            }
            else
            {
                check_result_unsat();
            }
            m_assignment.copy(m_temp_assignment);
        }

        /**
         * @brief Check whether the final result is correct when find sat status
         * 1. First check if all constraints are satisfied
         * 2. Second check if sum is the best objective sum
         */
        void check_result_sat()
        {
            m_assignment.copy(m_best_obj_assignment);
            for (int i = 0; i < m_constraints.size(); i++)
            {
                long double sum = m_constraints[i]->m_poly.calculate_sum(m_assignment);
                auto res = check_constraint_sat_using_sum(m_constraints[i], sum);
                if (!(res.first && res.second))
                {
                    UNREACHABLE();
                }
            }
            long double temp_sum = 0;
            for (int i = 0; i < m_vars.size(); i++)
            {
                temp_sum += m_assignment.value(i) * m_vars[i].m_coeff;
            }
            ASSERT(temp_sum == m_best_object_sum);
        }

        /**
         * @brief Check whether the final result is correct when not find sat status (best unsat)
         * check if unsat num is equal to unsat constraints' num
         */
        void check_result_unsat()
        {
            m_assignment.copy(m_best_unsat_assignment);
            int unsat_num = 0;
            for (int i = 0; i < m_constraints.size(); i++)
            {
                long double sum = m_constraints[i]->m_poly.calculate_sum(m_assignment);
                auto res = check_constraint_sat_using_sum(m_constraints[i], sum);
                if (!(res.first && res.second))
                {
                    unsat_num++;
                }
            }
            ASSERT(unsat_num == m_best_unsat_num);
        }

        /**
         * @brief Generate initialize assignment
         * 1. enable greedy initialize
         * 2. basic initialize: all var's are set to zero
         */
        void initialize_assignment()
        {
            if (enable_initialize_greedy)
            {
                initialize_assignment_greedy();
            }
            else
            {
                initialize_assignment_ordinary();
            }
        }

        void initialize_assignment_ordinary()
        {
            for (int i = 0; i < m_vars.size(); i++)
            {
                m_assignment.set(i, 0);
            }
        }

        void store_var_cast()
        {
            for (auto idx : m_cast_cons)
            {
                auto cons = m_constraints[idx];
                cons->m_poly.store_vars_and_generate_factor();
                for (int v : cons->m_poly.m_vars)
                {
                    m_vars[v].m_cast_id = idx;
                }
                cast_cons_left[idx] = cons->m_upper_bound;
            }
        }

        void initialize_assignment_greedy()
        {
            store_var_cast();
            for (auto idx : m_demand_cons)
            {
                auto cons = m_constraints[idx];
                cons->m_poly.store_vars_and_generate_factor();
                ASSERT(cons->m_lower);
                long double curr_sum = 0;
                long double demand_amount = cons->m_lower_bound;
                for (int v : cons->m_poly.m_vars)
                {
                    auto cast_id = m_vars[v].m_cast_id;
                    long double left_cast = cast_cons_left[cast_id];
                    m_assignment.set(v, std::min(m_vars[v].m_ub, std::min(demand_amount - curr_sum, left_cast)));
                    curr_sum += m_assignment.value(v);
                    cast_cons_left[cast_id] -= m_assignment.value(v);
                }
            }
        }

        static std::string bool2str(bool b)
        {
            return b ? "true" : "false";
        }

        static std::string boolpair2str(bool_pair b)
        {
            return bool2str(b.first) + ", " + bool2str(b.second);
        }

        /**
         * @brief Initialize constraints' status
         * 1. constraint's polynomial sum
         * 2. lower sat and upper sat
         * 3. lower bound and upper bound
         * 4. overall unsat and bound status
         */
        void initialize_constraints()
        {
            m_unsat_constraints.clear();
            for (auto cons : m_constraints)
            {
                long double curr_sum = cons->m_poly.calculate_sum(m_assignment);
                cons->m_poly.set_sum(curr_sum);
                bool_pair curr_sat = check_constraint_sat_using_sum(cons, curr_sum);
                bool overall_sat = curr_sat.first && curr_sat.second;
                cons->set_lower_sat(curr_sat.first);
                cons->set_upper_sat(curr_sat.second);
                cons->set_sat_status(overall_sat);
                // std::cout << "name: " << cons->m_name << std::endl;
                // std::cout << "sum: " << curr_sum << ", lb: " << cons->m_lower_bound << ", ub: " << cons->m_upper_bound << std::endl;
                // std::cout << boolpair2str(curr_sat) << std::endl;
                bool curr_bound = check_constraint_bound_using_sum(cons, curr_sum);
                cons->set_bound_status(curr_bound);
                if (!overall_sat)
                {
                    m_unsat_constraints.insert(cons->get_index());
                }
                if (!curr_bound)
                {
                    m_unbounded_constraints.insert(cons->get_index());
                }
                cons->m_poly.store_vars_and_generate_factor();
                cons->m_poly.generate_factor_values(m_assignment);
                cons->m_poly.generate_var_neighbours();
            }
            is_sat_status = m_unsat_constraints.empty();
        }

        void initialize_var_constraints()
        {
            for (auto cons : m_constraints)
            {
                for (int v : cons->m_poly.m_vars)
                {
                    m_vars[v].add_constraint(cons->get_index());
                }
            }
        }

        /**
         * @brief Initialize objective sum
         *
         */
        void initialize_objective()
        {
            // Generate factor and intercept vars for constraints containing objective vars
            // for(auto cons: m_constraints) {
            //     cons->m_poly.store_factor_and_intercept_all_vars();
            // }
            // Calculate objective sum
            m_object_sum = 0;
            for (auto v : m_vars)
            {
                if (v.m_coeff != 0)
                {
                    m_object_sum += v.m_coeff * m_assignment.value(v.m_index);
                }
            }
            m_object_weight = 0;
        }

        bool_pair check_constraint_sat_using_sum(constraint const *cons, long double sum) const
        {
            bool lower_sat = !cons->m_lower || sum >= cons->m_lower_bound;
            bool upper_sat = !cons->m_upper || sum <= cons->m_upper_bound;
            return std::make_pair(lower_sat, upper_sat);
        }

        /**
         * @brief Given a constraint and its sum, check whether the constraint is bounded
         *
         */
        bool check_constraint_bound_using_sum(constraint const *cons, long double sum) const
        {
            if (cons->m_lower && cons->m_upper)
            {
                return cons->m_lower_bound == cons->m_upper_bound && sum == cons->m_lower_bound;
            }
            else if (cons->m_lower && !cons->m_upper)
            {
                return cons->m_lower_bound == sum;
            }
            else if (!cons->m_lower && cons->m_upper)
            {
                return cons->m_upper_bound == sum;
            }
            else
            {
                UNREACHABLE();
                return false;
            }
        }

        /**
         * @brief Insert operations from a single constraint
         *
         */
        void insert_operation_from_single_unbound_constraint(constraint *cons)
        {
            if (cons->m_lower)
            {
                insert_operation_from_single_constraint_lower(cons);
            }
            if (cons->m_upper)
            {
                insert_operation_from_single_constraint_upper(cons);
            }
        }

        void insert_operation_from_single_constraint_lower(constraint *cons)
        {
            long double sum = cons->m_poly.get_sum();
            /*
                p * x + q == sum, p * (x + shift) + q == lb
                -------------------------------------------
                shift == (lb - sum) / p
            */
            long double diff = cons->m_lower_bound - sum;
            for (int v : cons->m_poly.m_vars)
            {
                long double p_val = cons->m_poly.m_var_factor_value.at(v);
                if (p_val == 0)
                {
                    continue;
                }
                long double real_shift = diff / p_val;
                long double shift1 = std::ceil(real_shift), shift2 = std::floor(real_shift);
                long double sum1 = sum + p_val * shift1, sum2 = sum + p_val * shift2;
                bool_pair sat1 = check_constraint_sat_using_sum(cons, sum1);
                long double shift = sat1.first ? shift1 : shift2;
                if (check_var_shift(v, shift))
                {
                    m_operation_table.insert(v, shift);
                }
            }
        }

        void insert_operation_from_single_constraint_upper(constraint *cons)
        {
            long double sum = cons->m_poly.get_sum();
            /*
                p * x + q == sum, p * (x + shift) + q == ub
                -------------------------------------------
                shift == (ub - sum) / p
            */
            long double diff = cons->m_upper_bound - sum;
            for (int v : cons->m_poly.m_vars)
            {
                long double p_val = cons->m_poly.m_var_factor_value.at(v);
                if (p_val == 0)
                {
                    continue;
                }
                long double real_shift = diff / p_val;
                long double shift1 = std::ceil(real_shift), shift2 = std::floor(real_shift);
                long double sum1 = sum + p_val * shift1, sum2 = sum + p_val * shift2;
                bool_pair sat1 = check_constraint_sat_using_sum(cons, sum1);
                long double shift = sat1.second ? shift1 : shift2;
                if (check_var_shift(v, shift))
                {
                    m_operation_table.insert(v, shift);
                }
            }
        }

        /**
         * @brief Given a var index and its change value, check whtether the change value is legal
         * 1. ignore shift == 0
         * 2. consider tabu
         * 3. post value still in var's bound
         */
        bool check_var_shift(int idx, long double shift) const
        {
            if (!is_number(shift) || !is_finite_number(shift))
            {
                return false;
            }
            if (shift == 0)
                return false;
            if (enable_tabu)
            {
                if ((shift > 0 && m_step <= m_vars[idx].m_pos_step) || (shift < 0 && m_step <= m_vars[idx].m_neg_step))
                {
                    return false;
                }
            }
            // return true;
            // if(enable_bound_constraints) {
            //     return true;
            // }
            long double post_val = m_assignment.value(idx) + shift;
            bool lower_sat = !m_vars[idx].m_lower || post_val >= m_vars[idx].m_lb;
            bool upper_sat = !m_vars[idx].m_upper || post_val <= m_vars[idx].m_ub;
            return lower_sat && upper_sat;
            // return lower_sat;
        }

        /**
         * @brief Insert operation from unbound constraints which contain at least one objective var
         *
         */
        void insert_operation_from_unbound_constraints_objective()
        {
            m_operation_table.clear();
            int_table m_cons_table;
            // Loop constraints linked by vars in objective
            for (auto v : m_vars_in_obj)
            {
                for (auto idx : m_vars[v].m_constraints)
                {
                    if (m_unbounded_constraints.count(idx) == 0 || m_cons_table.count(idx) != 0)
                    {
                        continue;
                    }
                    auto cons = m_constraints[idx];
                    m_cons_table.insert(idx);
                    insert_operation_from_single_unbound_constraint(cons);
                }
            }
        }

        /**
         * @brief Insert operation from a single unsat constraint
         *
         */
        void insert_operation_from_single_unsat_constraint(constraint *cons)
        {
            long double sum = cons->m_poly.get_sum();
            if (cons->m_lower && !cons->m_upper)
            {
                ASSERT(sum < cons->m_lower_bound);
                insert_operation_from_single_constraint_lower(cons);
            }
            else if (!cons->m_lower && cons->m_upper)
            {
                ASSERT(sum > cons->m_upper_bound);
                insert_operation_from_single_constraint_upper(cons);
            }
            else if (cons->m_lower && cons->m_upper)
            {
                if (sum < cons->m_lower_bound)
                {
                    insert_operation_from_single_constraint_lower(cons);
                }
                else if (sum > cons->m_upper_bound)
                {
                    insert_operation_from_single_constraint_upper(cons);
                }
                else
                {
                    UNREACHABLE();
                }
            }
            else
            {
                UNREACHABLE();
            }
        }

        /**
         * @brief Insert operations from all unsat constraints this step
         *
         */
        void insert_operation_from_unsat_constraints()
        {
            m_operation_table.clear();
            for (auto idx : m_unsat_constraints)
            {
                constraint *cons = m_constraints[idx];
                insert_operation_from_single_unsat_constraint(cons);
            }
        }

        void select_best_operation(int &index, long double &shift, long double &score, bool is_unsat = false)
        {
            int cnt;
            bool BMS;
            // std::cout << "m_operation_table.size: " << m_operation_table.size() << std::endl;
            if (m_operation_table.size() <= BMS_VAL)
            {
                cnt = m_operation_table.size();
                BMS = false;
            }
            else
            {
                cnt = BMS_VAL;
                BMS = true;
            }
            int curr_var;
            long double curr_shift, curr_score;
            score = INT32_MIN;
            int size = m_operation_table.size();
            for (int i = 0; i < cnt; i++)
            {
                if (BMS)
                {
                    int rand_index = rand() % (size - i);
                    curr_var = m_operation_table.m_vars[rand_index];
                    curr_shift = m_operation_table.m_shifts[rand_index];
                    m_operation_table.m_vars[rand_index] = m_operation_table.m_vars[size - i - 1];
                    m_operation_table.m_shifts[rand_index] = m_operation_table.m_shifts[size - i - 1];
                }
                else
                {
                    curr_var = m_operation_table.m_vars[i];
                    curr_shift = m_operation_table.m_shifts[i];
                }
                // if(is_unsat && std::find(m_vars_in_obj.begin(),m_vars_in_obj.end(),curr_var)!=m_vars_in_obj.end())
                // {
                //     continue;
                // }
                curr_score = calculate_move_score(curr_var, curr_shift);
                if (curr_score > score)
                {
                    score = curr_score;
                    index = curr_var;
                    shift = curr_shift;
                }
            }
        }

        /**
         * @brief move score = hard score + soft score
         */
        long double calculate_move_score(int m_var, long double m_shift)
        {
            return hard_score(m_var, m_shift) + soft_score(m_var, m_shift);
        }

        /**
         * @brief unsat clause's overall weight change value
         */
        long double hard_score(int m_var, long double m_shift)
        {
            long double score = 0;
            for (auto cons : m_constraints)
            {
                if (cons->m_poly.m_vars.count(m_var) == 0)
                {
                    continue;
                }
                bool pre_sat = cons->get_sat_status();
                long double pre_sum = cons->m_poly.get_sum();
                long double post_sum = pre_sum + cons->m_poly.m_var_factor_value.at(m_var) * m_shift;
                bool_pair post_sat_pair = check_constraint_sat_using_sum(cons, post_sum);
                bool post_sat = post_sat_pair.first && post_sat_pair.second;
                // case I. sat --> unsat
                if (pre_sat && !post_sat)
                {
                    if (cons->m_lower && cons->m_upper)
                    {
                        if (post_sum < cons->m_lower_bound)
                        {
                            score -= cons->get_lower_weight() * 1.0;
                        }
                        else if (post_sum > cons->m_upper_bound)
                        {
                            score -= cons->get_upper_weight() * 1.0;
                        }
                        else
                        {
                            UNREACHABLE();
                        }
                    }
                    else if (!cons->m_lower && cons->m_upper)
                    {
                        ASSERT(post_sum > cons->m_upper_bound);
                        score -= cons->get_upper_weight() * 1.0;
                    }
                    else if (cons->m_lower && !cons->m_upper)
                    {
                        ASSERT_CODE(post_sum < cons->m_lower_bound,
                                    std::cout << "post sum: " << post_sum << std::endl;
                                    std::cout << "lower bound: " << cons->m_lower_bound << std::endl;
                                    std::cout << "pre sum: " << pre_sum << std::endl;
                                    std::cout << "coeff: " << cons->m_poly.m_var_factor_value.at(m_var) << std::endl;
                                    std::cout << "shift: " << m_shift << std::endl;);
                        score -= cons->get_lower_weight() * 1.0;
                    }
                    else
                    {
                        UNREACHABLE();
                    }
                }
                else if (!pre_sat && post_sat)
                { // case II. unsat --> sat
                    if (cons->m_lower && cons->m_upper)
                    {
                        if (pre_sum < cons->m_lower_bound)
                        {
                            score += cons->get_lower_weight() * 1.0;
                        }
                        else if (pre_sum > cons->m_upper_bound)
                        {
                            score += cons->get_upper_weight() * 1.0;
                        }
                        else
                        {
                            UNREACHABLE();
                        }
                    }
                    else if (!cons->m_lower && cons->m_upper)
                    {
                        ASSERT(pre_sum > cons->m_upper_bound);
                        score += cons->get_upper_weight() * 1.0;
                    }
                    else if (cons->m_lower && !cons->m_upper)
                    {
                        ASSERT(pre_sum < cons->m_lower_bound);
                        score += cons->get_lower_weight() * 1.0;
                    }
                    else
                    {
                        UNREACHABLE();
                    }
                }
                else if (!pre_sat && !post_sat)
                { // case III. unsat --> unsat, only consider two bounds' case
                    if (cons->m_lower && cons->m_upper)
                    {
                        if (pre_sum < cons->m_lower_bound && post_sum > cons->m_upper_bound)
                        { // lower unsat --> upper unsat
                            score += cons->get_lower_weight() * 1.0;
                            score -= cons->get_upper_weight() * 1.0;
                        }
                        else if (pre_sum > cons->m_upper_bound && post_sum < cons->m_lower_bound)
                        { // upper unsat --> lower unsat
                            score += cons->get_upper_weight() * 1.0;
                            score -= cons->get_lower_weight() * 1.0;
                        }
                        else
                        { // other cases, don't consider
                        }
                    }
                }
                else
                { // case IV. sat --> sat, don't consider
                }
            } // finish loop
            return score;
        }

        /**
         * @brief Currently direction score
         *
         */
        long double soft_score(int m_var, long double m_shift) const
        {
            long double obj_move = m_vars[m_var].m_coeff * m_shift;
            return obj_move == 0 ? 0 : (obj_move < 0 ? 1.0 * m_object_weight : -1.0 * m_object_weight);
        }

        /**
         * @brief adjust var's factor value in all constraints which contain this var
         *
         */
        void adjust_factor_value(int v, constraint *cons)
        {
            if (cons->m_poly.contains(v))
            {
                cons->m_poly.generate_factor_values_core(v, m_assignment);
            }
        }

        /**
         * @brief move a var with change value
         *
         */
        void execute_critical_move(int m_var, long double m_shift)
        {
            TRACE(tout << "execute critical move\n";
                  tout << m_var << ", " << m_vars[m_var].m_name << " -> ";
                  display_sign_number(tout, m_shift) << " (" << m_assignment.value(m_var) + m_shift << ") " << std::endl;);
            std::cout << "  var:  " << m_var << "   m_shift  " << m_shift << std::endl;
            if (m_shift == 0)
            {
                return;
            }
            for (auto idx : m_vars[m_var].m_constraints)
            {
                auto cons = m_constraints[idx];
                bool pre_sat = cons->get_sat_status(), pre_bound = cons->get_bound_status();
                long double sum = cons->m_poly.get_sum() + cons->m_poly.m_var_factor_value.at(m_var) * m_shift;
                bool_pair post_sat_pair = check_constraint_sat_using_sum(cons, sum);
                bool post_sat = post_sat_pair.first && post_sat_pair.second;
                bool post_bound = check_constraint_bound_using_sum(cons, sum);
                if (pre_sat && !post_sat)
                { // sat --> unsat
                    m_unsat_constraints.insert(cons->get_index());
                }
                else if (!pre_sat && post_sat)
                { // unsat --> sat
                    m_unsat_constraints.erase(cons->get_index());
                }
                if (pre_bound && !post_bound)
                { // bound --> unbound
                    m_unbounded_constraints.insert(cons->get_index());
                }
                else if (!pre_bound && post_bound)
                { // unbound --> bound
                    m_unbounded_constraints.erase(cons->get_index());
                }
                cons->set_sat_status(post_sat);
                cons->set_lower_sat(post_sat_pair.first);
                cons->set_upper_sat(post_sat_pair.second);
                cons->set_bound_status(post_bound);
                cons->m_poly.set_sum(sum);
            }
            m_assignment.shift(m_var, m_shift);
            if (enable_tabu)
            {
                if (m_shift > 0)
                {
                    m_vars[m_var].m_neg_step = m_step + tabu_const + std::rand() % tabu_rand;
                }
                else if (m_shift < 0)
                {
                    m_vars[m_var].m_pos_step = m_step + tabu_const + std::rand() % tabu_rand;
                }
            }
            is_sat_status = m_unsat_constraints.empty();
            m_object_sum += m_vars[m_var].m_coeff * m_shift;
            for (auto idx : m_vars[m_var].m_constraints)
            {
                auto cons = m_constraints[idx];
                for (int nv : cons->m_poly.m_var_neighbours[m_var])
                {
                    adjust_factor_value(nv, cons);
                }
            }
            TRACE(
                display_assignment(tout, m_assignment);
                tout << "object sum: " << m_object_sum << std::endl;);
        }

        /**
         * @brief random select num elements from st into vec
         *
         */
        void random_select_from_set(int_table const &st, int num, int_vector &vec) const
        {
            vec.clear();
            int_table random_indices;
            for (int i = 0; i < num; i++)
            {
                random_indices.insert(rand() % st.size());
            }
            int cnt = 0;
            for (auto it = st.begin(); it != st.end(); it++)
            {
                if (random_indices.count(cnt) != 0)
                {
                    vec.push_back(*it);
                }
                cnt++;
            }
        }

        /*
            TODO: move var to lb or ub according to constraint's move direction
        */
        void random_walk_sat()
        {
            m_operation_table.clear();
            if (m_unbounded_constraints.empty())
            {
                return;
            }
            int_vector cons_indices;
            random_select_from_set(m_unbounded_constraints, random_walk_num, cons_indices);
            for (auto i : cons_indices)
            {
                auto cons = m_constraints[i];
                insert_operation_from_single_unbound_constraint(cons);
            }
            if (m_operation_table.empty())
            { // no operation
                no_operation_random_walk_unbound();
            }
            else
            { // choose the best operation
                int v_idx;
                long double v_shift, score;
                select_best_operation(v_idx, v_shift, score);
                execute_critical_move(v_idx, v_shift);
            }
        }

        void random_walk_unsat()
        {
            m_operation_table.clear();
            ASSERT(!m_unsat_constraints.empty());
            int_vector cons_indices;
            random_select_from_set(m_unsat_constraints, random_walk_num, cons_indices);
            for (auto i : cons_indices)
            {
                auto cons = m_constraints[i];
                insert_operation_from_single_unsat_constraint(cons);
            }
            if (m_operation_table.empty())
            { // no operation
                no_operation_random_walk_unsat();
            }
            else
            { // choose the best operation
                int v_idx;
                long double v_shift, score;
                select_best_operation(v_idx, v_shift, score);
                execute_critical_move(v_idx, v_shift);
            }
        }

        void no_operation_random_walk_unbound()
        {
            /*
                Random select an unbounded constraint, random select on var
                Move to its upper bound or lower bound
            */
            int_vector r_vec;
            random_select_from_set(m_unbounded_constraints, 1, r_vec);
            ASSERT(r_vec.size() == 1);
            auto cons = m_constraints[r_vec[0]];
            random_select_from_set(cons->m_poly.m_vars, 1, r_vec);
            ASSERT(r_vec.size() == 1);
            int v_idx = r_vec[0];
            long double shift;
            if (m_vars[v_idx].m_lower && m_vars[v_idx].m_upper)
            {
                int r = rand() % 2;
                if (r == 0)
                {
                    shift = m_vars[v_idx].m_lb - m_assignment.value(v_idx);
                }
                else
                {
                    shift = m_vars[v_idx].m_ub - m_assignment.value(v_idx);
                }
            }
            else if (!m_vars[v_idx].m_lower && m_vars[v_idx].m_upper)
            {
                shift = m_vars[v_idx].m_ub - m_assignment.value(v_idx);
            }
            else if (m_vars[v_idx].m_lower && !m_vars[v_idx].m_upper)
            {
                shift = m_vars[v_idx].m_lb - m_assignment.value(v_idx);
            }
            else
            {
                return;
            }
            execute_critical_move(v_idx, shift);
        }

        void no_operation_random_walk_unsat()
        {
            /*
                Random select an unsat constraint, random select on var
                Move to its upper bound or lower bound
            */
            int_vector r_vec;
            random_select_from_set(m_unsat_constraints, 1, r_vec);
            ASSERT(r_vec.size() == 1);
            auto cons = m_constraints[r_vec[0]];
            random_select_from_set(cons->m_poly.m_vars, 1, r_vec);
            ASSERT(r_vec.size() == 1);
            int v_idx = r_vec[0];
            long double shift;
            if (m_vars[v_idx].m_lower && m_vars[v_idx].m_upper)
            {
                int r = rand() % 2;
                if (r == 0)
                {
                    shift = m_vars[v_idx].m_lb - m_assignment.value(v_idx);
                }
                else
                {
                    shift = m_vars[v_idx].m_ub - m_assignment.value(v_idx);
                }
            }
            else if (!m_vars[v_idx].m_lower && m_vars[v_idx].m_upper)
            {
                shift = m_vars[v_idx].m_ub - m_assignment.value(v_idx);
            }
            else if (m_vars[v_idx].m_lower && !m_vars[v_idx].m_upper)
            {
                shift = m_vars[v_idx].m_lb - m_assignment.value(v_idx);
            }
            else
            {
                return;
            }
            execute_critical_move(v_idx, shift);
        }

        void update_weight()
        {
            for (auto idx : m_unsat_constraints)
            {
                auto cons = m_constraints[idx];
                long double sum = cons->m_poly.get_sum();
                if (cons->m_lower && cons->m_upper)
                {
                    if (sum < cons->m_lower_bound)
                    {
                        cons->inc_lower_weight();
                        m_sum_weight++;
                    }
                    else if (sum > cons->m_upper_bound)
                    {
                        cons->inc_upper_weight();
                        m_sum_weight++;
                    }
                    else
                    {
                        UNREACHABLE();
                    }
                }
                else if (!cons->m_lower && cons->m_upper)
                {
                    ASSERT(sum > cons->m_upper_bound);
                    cons->inc_upper_weight();
                    m_sum_weight++;
                }
                else if (cons->m_lower && !cons->m_upper)
                {
                    ASSERT(sum < cons->m_lower_bound);
                    cons->inc_lower_weight();
                    m_sum_weight++;
                }
                else
                {
                    UNREACHABLE();
                }
            }
            m_avg_weight = m_sum_weight * 1.0 / m_constraints.size();
            /**
             * Two Prerequisites for updating objecive weight
             * 1. obj_weight < avg_weight + GAP
             * 2. not updating best objective value
             */
            // if(m_object_weight < m_avg_weight + weight_gap && find_sat_status && m_object_sum >= m_best_object_sum) {
            if (m_object_weight < 100 && find_sat_status && m_object_sum >= m_best_object_sum)
            {
                m_object_weight++;
            }
        }

        void update_best_information()
        {
            // update unsat information
            if (m_unsat_constraints.size() < m_best_unsat_num)
            {
                m_best_unsat_num = m_unsat_constraints.size();
                m_best_unsat_assignment.copy(m_assignment);
            }
            // update best objective information
            if (is_sat_status)
            {
                if (!find_sat_status)
                {
                    m_object_weight = 1;
                    TRACE(
                        tout << "find feasible assignment\n";);
                }
                if (!find_sat_status || (find_sat_status && m_object_sum < m_best_object_sum))
                {
                    m_best_object_sum = m_object_sum;
                    TRACE(display_assignment(tout, m_assignment););
                    m_best_obj_assignment.copy(m_assignment);
                }
                find_sat_status = true;
            }
        }

        std::ostream &display_all_constraints(std::ostream &out) const
        {
            out << "display all constraints\n";
            for (auto cons : m_constraints)
            {
                display_constraint(out, *cons);
            }
            out << std::endl;
            return out;
        }

        std::ostream &display_all_constraints_with_weight(std::ostream &out) const
        {
            out << "display all constraints\n";
            for (auto cons : m_constraints)
            {
                display_constraint_with_weight(out, *cons);
            }
            out << std::endl;
            return out;
        }

        std::ostream &display_sign_number(std::ostream &out, long double x) const
        {
            if (x > 0)
            {
                out << "+";
            }
            out << x;
            return out;
        }

        std::ostream &display_operation_table(std::ostream &out, operation_table const &table) const
        {
            out << "display operation table\n";
            for (int i = 0; i < table.m_vars.size(); i++)
            {
                out << m_vars[table.m_vars[i]].m_name << " -> ";
                display_sign_number(out, table.m_shifts[i]) << " (" << m_assignment.value(table.m_vars[i]) + table.m_shifts[i] << ") " << std::endl;
            }
            out << std::endl;
            return out;
        }

        std::ostream &display_best_choice(std::ostream &out, int idx, long double shift, long double score) const
        {
            out << "display best choice\n";
            if (idx == INT32_MAX)
            {
                out << "no information\n";
                out << std::endl;
                return out;
            }
            out << idx << ", " << m_vars[idx].m_name << " -> ";
            display_sign_number(out, shift) << " (" << m_assignment.value(idx) + shift << ") " << ", score: " << score << std::endl;
            out << std::endl;
            return out;
        }

        std::ostream &display_assignment(std::ostream &out, assignment const &ass) const
        {
            if (ass.empty())
            {
                out << "empty assignment\n";
            }
            else
            {
                // std::cout <<"var: size"<< m_vars.size() << std::endl;
                for (int i = 0; i < m_vars.size(); i++)
                {
                    out << m_vars[i].m_name << " -> " << ass.value(i) << std::endl;
                }
            }
            return out;
        }

        std::ostream &display_unsat_constraints(std::ostream &out) const
        {
            out << "display unsat constraints\n";
            out << "#unsat: " << m_unsat_constraints.size() << std::endl;
            for (auto idx : m_unsat_constraints)
            {
                constraint const *cons = m_constraints[idx];
                if (!cons->get_lower_sat())
                {
                    display_constraint_lower(out, *cons);
                }
                if (!cons->get_upper_sat())
                {
                    display_constraint_upper(out, *cons);
                }
            }
            out << std::endl;
            return out;
        }

        std::ostream &display_find_sat(std::ostream &out) const
        {
            if (find_sat_status)
            {
                out << "find feasible assignment\n";
                display_assignment(out, m_best_obj_assignment);
                out << "best object sum: " << m_best_object_sum << std::endl;
            }
            else
            {
                out << "no feasible assignment\n";
                display_assignment(out, m_best_unsat_assignment);
                out << "best unsat num: " << m_best_unsat_num << std::endl;
            }
            return out;
        }

        void insert_and_select_lift_moves(int &v_idx, long double &shift, long double &score)
        {
            m_operation_table.clear();
            long double curr_shift;
            long double curr_score;
            score = INT32_MIN;
            for (int curr_var : m_vars_in_obj)
            {
                insert_lift_move(curr_var, curr_shift, curr_score);
                if (enable_tabu)
                {
                    if ((curr_shift > 0 && m_step <= m_vars[curr_var].m_pos_step) || (curr_shift < 0 && m_step <= m_vars[curr_var].m_neg_step))
                    {
                        continue;
                    }
                }
                if (curr_score > score)
                {
                    score = curr_score;
                    v_idx = curr_var;
                    shift = curr_shift;
                }
            }
        }

        void insert_lift_move(int m_var, long double &shift, long double &score)
        {
            ASSERT(m_vars[m_var].m_coeff != 0);
            m_vars[m_var].m_crucial_constraint.clear();
            int val;
            bool inf, exist;
            if (m_vars[m_var].m_coeff > 0)
            {
                get_smallest_feasible_int_value(m_var, val, inf, exist);
            }
            else
            {
                get_largest_feasible_int_value(m_var, val, inf, exist);
            }
            ASSERT(exist);
            if (inf)
            {
                UNREACHABLE();
            }
            else
            {
                shift = val - m_assignment.value(m_var);
                score = -1.0 * m_vars[m_var].m_coeff * shift;
                ASSERT(is_in_bound(m_var, val));
            }
        }

        inline bool is_in_bound(int m_var, double m_val) const
        {
            return (!m_vars[m_var].m_lower || m_val >= m_vars[m_var].m_lb) && (!m_vars[m_var].m_upper || m_val <= m_vars[m_var].m_ub);
        }

        void insert_and_select_vnd_lift_moves(int &v1, int &v2, long double &shift1, long double &shift2, long double &score)
        {
            score = INT32_MIN;
            long double curr_shift1, curr_shift2, curr_score;
            // Loop all var-pairs in objective
            for (int curr_v1 : m_vars_in_obj)
            {
                for (int curr_v2 : m_vars_in_obj)
                {
                    if (curr_v1 < curr_v2)
                    {
                        // 0 means the distance from cur assignment upper bound, 1 means upper bound, 2 means 1：1, 3 means random
                        if (insert_vnd_lift_move(curr_v1, curr_v2, curr_shift1, curr_shift2, curr_score, 0))
                        {
                            if (curr_score > score)
                            {
                                score = curr_score;
                                v1 = curr_v1;
                                v2 = curr_v2;
                                shift1 = curr_shift1;
                                shift2 = curr_shift2;
                            }
                        }
                        if (insert_vnd_lift_move(curr_v1, curr_v2, curr_shift1, curr_shift2, curr_score, 1))
                        {
                            if (curr_score > score)
                            {
                                score = curr_score;
                                v1 = curr_v1;
                                v2 = curr_v2;
                                shift1 = curr_shift1;
                                shift2 = curr_shift2;
                            }
                        }
                        if (insert_vnd_lift_move(curr_v1, curr_v2, curr_shift1, curr_shift2, curr_score, 2))
                        {
                            if (curr_score > score)
                            {
                                score = curr_score;
                                v1 = curr_v1;
                                v2 = curr_v2;
                                shift1 = curr_shift1;
                                shift2 = curr_shift2;
                            }
                        }
                        // if(insert_vnd_lift_move(curr_v1, curr_v2, curr_shift1, curr_shift2, curr_score, 3)) {
                        //     if(curr_score > score) {
                        //         score = curr_score;
                        //         v1 = curr_v1; v2 = curr_v2;
                        //         shift1 = curr_shift1; shift2 = curr_shift2;
                        //     }
                        // }
                        if (insert_vnd_lift_move(curr_v1, curr_v2, curr_shift1, curr_shift2, curr_score, 4))
                        {
                            if (curr_score > score)
                            {
                                score = curr_score;
                                v1 = curr_v1;
                                v2 = curr_v2;
                                shift1 = curr_shift1;
                                shift2 = curr_shift2;
                            }
                        }
                    }
                }
            }
        }

        /**
         * e.g. obj: x + 2y
         * then x = x + 2t, y = y + t
         */
        bool insert_vnd_lift_move(int v1, int v2, long double &shift1, long double &shift2, long double &score, int vnd_mode)
        {
            int_table m_looped_cons;
            auto coeff1 = m_vars[v1].m_coeff, coeff2 = m_vars[v2].m_coeff;
            int t_coeff1, t_coeff2;
            switch (vnd_mode)
            {
            case 0:
            {
                t_coeff1 = m_vars[v1].m_ub - m_assignment.value(v1);
                t_coeff2 = m_vars[v2].m_ub - m_assignment.value(v2);
                break;
            }
            case 1:
            {
                t_coeff1 = m_vars[v1].m_ub;
                t_coeff2 = m_vars[v2].m_ub;
                break;
            }
            case 2:
            {
                t_coeff1 = 1;
                t_coeff2 = 1;
                break;
            }
            case 3:
            {
                t_coeff1 = 1;
                t_coeff2 = 1;
                break;
            }
            case 4:
            {
                t_coeff1 = m_assignment.value(v1);
                t_coeff2 = m_assignment.value(v2);
                break;
            }
            default:
                std::cout << "error" << std::endl;
            }
            if (t_coeff1 == 0 || t_coeff2 == 0)
                return false;
            // coeff1 * v1 + coeff2 * v2
            // coeff1 * (v1 + t_coeff1 * t) + coeff2 * (v2 + t_coeff2 * t)
            auto t_coeff = coeff1 * t_coeff1 + coeff2 * t_coeff2;
            if (t_coeff == 0)
                return false;
            ASSERT(t_coeff != 0);
            // Calculate feasible set of t
            // Currently we do not consider var's bound because we have bound constraints
            bool inf = true;
            long double val = t_coeff > 0 ? INT32_MIN : INT32_MAX;
            for (int idx : m_vars[v1].m_constraints)
            {
                auto cons = m_constraints[idx];
                bool curr_exist, curr_inf;
                long double curr_val;
                if (t_coeff > 0)
                {
                    get_smallest_vnd_feasible_int_cons(cons, v1, v2, t_coeff1, t_coeff2, curr_exist, curr_inf, curr_val);
                    if (!curr_exist)
                    {
                        return false;
                    }
                    inf = inf && curr_inf;
                    if (!curr_inf)
                    {
                        val = std::max(val, curr_val);
                    }
                }
                else
                {
                    get_largest_vnd_feasible_int_cons(cons, v1, v2, t_coeff1, t_coeff2, curr_exist, curr_inf, curr_val);
                    if (!curr_exist)
                    {
                        return false;
                    }
                    inf = inf && curr_inf;
                    if (!curr_inf)
                    {
                        val = std::min(val, curr_val);
                    }
                }
                m_looped_cons.insert(idx);
            }
            for (int idx : m_vars[v2].m_constraints)
            {
                if (m_looped_cons.count(idx) != 0)
                {
                    continue;
                }
                auto cons = m_constraints[idx];
                bool curr_exist, curr_inf;
                long double curr_val;
                if (t_coeff > 0)
                {
                    get_smallest_vnd_feasible_int_cons(cons, v1, v2, t_coeff1, t_coeff2, curr_exist, curr_inf, curr_val);
                    if (!curr_exist)
                    {
                        return false;
                    }
                    inf = inf && curr_inf;
                    if (!curr_inf)
                    {
                        val = std::max(val, curr_val);
                    }
                }
                else
                {
                    get_largest_vnd_feasible_int_cons(cons, v1, v2, t_coeff1, t_coeff2, curr_exist, curr_inf, curr_val);
                    if (!curr_exist)
                    {
                        return false;
                    }
                    inf = inf && curr_inf;
                    if (!curr_inf)
                    {
                        val = std::min(val, curr_val);
                    }
                }
            }
            if (inf)
            {
                UNREACHABLE();
                return false;
            }
            else
            {
                shift1 = val * t_coeff1, shift2 = val * t_coeff2;
                score = -1.0 * t_coeff * val;
                return true;
            }
        }

        void get_smallest_vnd_feasible_int_cons(constraint *cons, int v1, int v2, long double t_coeff1, long double t_coeff2, bool &exist, bool &inf, long double &val)
        {
            if (cons->m_lower && cons->m_upper)
            {
                bool exist1, inf1, exist2, inf2;
                long double val1, val2;
                get_smallest_vnd_feasible_int_cons_lower(cons, v1, v2, t_coeff1, t_coeff2, exist1, inf1, val1);
                get_smallest_vnd_feasible_int_cons_upper(cons, v1, v2, t_coeff1, t_coeff2, exist2, inf2, val2);
                if (!(exist1 && exist2))
                {
                    exist = false;
                    return;
                }
                if (inf1 && inf2)
                {
                    inf = true;
                }
                else if (inf1 && !inf2)
                {
                    val = val2;
                }
                else if (!inf1 && inf2)
                {
                    val = val1;
                }
                else
                {
                    val = std::max(val1, val2);
                }
            }
            else if (!cons->m_lower && cons->m_upper)
            {
                get_smallest_vnd_feasible_int_cons_upper(cons, v1, v2, t_coeff1, t_coeff2, exist, inf, val);
            }
            else if (cons->m_lower && !cons->m_upper)
            {
                get_smallest_vnd_feasible_int_cons_lower(cons, v1, v2, t_coeff1, t_coeff2, exist, inf, val);
            }
            else
            {
                UNREACHABLE();
            }
        }

        // we ignore cases that v1 and v2 factor with each other
        void get_smallest_vnd_feasible_int_cons_lower(constraint *cons, int v1, int v2, long double t_coeff1, long double t_coeff2, bool &exist, bool &inf, long double &val)
        {
            auto sum = cons->m_poly.get_sum();
            auto lb = cons->m_lower_bound;
            auto diff = lb - sum;
            auto coeff1 = cons->m_poly.m_var_factor_value[v1], coeff2 = cons->m_poly.m_var_factor_value[v2];
            if (cons->m_poly.m_vars.count(v1) != 0 && cons->m_poly.m_vars.count(v2) != 0)
            {
                cons->m_poly.store_factor_and_intercept_vars(v1);
                if (cons->m_poly.m_var_vars[v1].first.count(v2) != 0)
                { // Factor with each other
                    exist = false;
                    return;
                }
                // coeff1 * v1 + coeff2 * v2 + other == sum
                // coeff1 * (v1 + t_coeff1 * t) + coeff2 * (v2 + t_coeff2 * t) + other >= lb
                // (coeff1 * t_coeff1 + coeff2 * t_coeff2) * t + sum >= lb
                // (coeff1 * t_coeff1 + coeff2 * t_coeff2) * t >= lb - sum
                // t >=< (lb - sum) / t_coeff
                auto t_coeff = coeff1 * t_coeff1 + coeff2 * t_coeff2;
                if (t_coeff == 0)
                {
                    if (sum >= lb)
                    {
                        exist = true;
                        inf = true;
                    }
                    else
                    {
                        exist = false;
                    }
                }
                else if (t_coeff > 0)
                {
                    exist = true;
                    inf = false;
                    val = std::ceil(diff / t_coeff);
                }
                else
                {
                    exist = true;
                    inf = true;
                }
            }
            else if (cons->m_poly.m_vars.count(v1) != 0 && cons->m_poly.m_vars.count(v2) == 0)
            {
                // coeff1 * v1 + other == sum
                // coeff1 * (v1 + t_coeff1 * t) + other >= lb
                // coeff1 * t_coeff1 * t >= lb - sum
                // t >=< (lb - sum) / t_coeff
                auto t_coeff = coeff1 * t_coeff1;
                if (t_coeff == 0)
                {
                    if (sum >= lb)
                    {
                        exist = true;
                        inf = true;
                    }
                    else
                    {
                        exist = false;
                    }
                }
                else if (t_coeff > 0)
                {
                    exist = true;
                    inf = false;
                    val = std::ceil(diff / t_coeff);
                }
                else
                {
                    exist = true;
                    inf = true;
                }
            }
            else if (cons->m_poly.m_vars.count(v1) == 0 && cons->m_poly.m_vars.count(v2) != 0)
            {
                // coeff2 * v2 + other == sum
                // coeff2 * (v2 + t_coeff2 * t) + other >= lb
                // coeff2 * t_coeff2 * t >= lb - sum
                // t >=< (lb - sum) / t_coeff
                auto t_coeff = coeff2 * t_coeff2;
                if (t_coeff == 0)
                {
                    if (sum >= lb)
                    {
                        exist = true;
                        inf = true;
                    }
                    else
                    {
                        exist = false;
                    }
                }
                else if (t_coeff > 0)
                {
                    exist = true;
                    inf = true;
                    val = std::ceil(diff / t_coeff);
                }
                else
                {
                    exist = true;
                    inf = true;
                }
            }
            else
            {
                exist = true;
                inf = true;
            }
        }

        void get_smallest_vnd_feasible_int_cons_upper(constraint *cons, int v1, int v2, long double t_coeff1, long double t_coeff2, bool &exist, bool &inf, long double &val)
        {
            auto sum = cons->m_poly.get_sum();
            auto ub = cons->m_upper_bound;
            auto diff = ub - sum;
            auto coeff1 = cons->m_poly.m_var_factor_value[v1], coeff2 = cons->m_poly.m_var_factor_value[v2];
            if (cons->m_poly.m_vars.count(v1) != 0 && cons->m_poly.m_vars.count(v2) != 0)
            {
                cons->m_poly.store_factor_and_intercept_vars(v1);
                if (cons->m_poly.m_var_vars[v1].first.count(v2) != 0)
                { // Factor with each other
                    exist = false;
                    return;
                }
                // coeff1 * v1 + coeff2 * v2 + other == sum
                // coeff1 * (v1 + t_coeff1 * t) + coeff2 * (v2 + t_coeff2 * t) + other  <= ub
                // t_coeff * t + sum <= ub
                auto t_coeff = coeff1 * t_coeff1 + coeff2 * t_coeff2;
                if (t_coeff == 0)
                {
                    if (sum <= ub)
                    {
                        exist = true;
                        inf = true;
                    }
                    else
                    {
                        exist = false;
                    }
                }
                else if (t_coeff > 0)
                {
                    exist = true;
                    inf = true;
                }
                else
                {
                    exist = true;
                    inf = false;
                    val = std::ceil(diff / t_coeff);
                }
            }
            else if (cons->m_poly.m_vars.count(v1) != 0 && cons->m_poly.m_vars.count(v2) == 0)
            {
                // coeff1 * v1 + other == sum
                // coeff1 * (v1 + t_coeff1 * t) + other <= ub
                // t_coeff * t + sum <= ub
                auto t_coeff = coeff1 * t_coeff1;
                if (t_coeff1 == 0)
                {
                    if (sum <= ub)
                    {
                        exist = true;
                        inf = true;
                    }
                    else
                    {
                        exist = false;
                    }
                }
                else if (t_coeff1 > 0)
                {
                    exist = true;
                    inf = true;
                }
                else
                {
                    exist = true;
                    inf = false;
                    val = std::ceil(diff / t_coeff);
                }
            }
            else if (cons->m_poly.m_vars.count(v1) == 0 && cons->m_poly.m_vars.count(v2) != 0)
            {
                // coeff2 * v2 + other == sum
                // coeff2 * (v2 + t_coeff2 * t) + other <= ub
                // t_coeff * t + sum <= ub
                auto t_coeff = coeff2 * t_coeff2;
                if (t_coeff == 0)
                {
                    if (sum <= ub)
                    {
                        exist = true;
                        inf = true;
                    }
                    else
                    {
                        exist = false;
                    }
                }
                else if (t_coeff > 0)
                {
                    exist = true;
                    inf = true;
                }
                else
                {
                    exist = true;
                    inf = false;
                    val = std::ceil(diff / t_coeff);
                }
            }
            else
            {
                exist = true;
                inf = true;
            }
        }

        void get_largest_vnd_feasible_int_cons(constraint *cons, int v1, int v2, long double t_coeff1, long double t_coeff2, bool &exist, bool &inf, long double &val)
        {
            if (cons->m_lower && cons->m_upper)
            {
                bool exist1, inf1, exist2, inf2;
                long double val1, val2;
                get_largest_vnd_feasible_int_cons_lower(cons, v1, v2, t_coeff1, t_coeff2, exist1, inf1, val1);
                get_largest_vnd_feasible_int_cons_upper(cons, v1, v2, t_coeff1, t_coeff2, exist2, inf2, val2);
                if (!(exist1 && exist2))
                {
                    exist = false;
                    return;
                }
                if (inf1 && inf2)
                {
                    inf = true;
                }
                else if (inf1 && !inf2)
                {
                    val = val2;
                }
                else if (!inf1 && inf2)
                {
                    val = val1;
                }
                else
                {
                    val = std::min(val1, val2);
                }
            }
            else if (!cons->m_lower && cons->m_upper)
            {
                get_largest_vnd_feasible_int_cons_upper(cons, v1, v2, t_coeff1, t_coeff2, exist, inf, val);
            }
            else if (cons->m_lower && !cons->m_upper)
            {
                get_largest_vnd_feasible_int_cons_lower(cons, v1, v2, t_coeff1, t_coeff2, exist, inf, val);
            }
            else
            {
                UNREACHABLE();
            }
        }

        void get_largest_vnd_feasible_int_cons_lower(constraint *cons, int v1, int v2, long double t_coeff1, long double t_coeff2, bool &exist, bool &inf, long double &val)
        {
            auto sum = cons->m_poly.get_sum();
            auto lb = cons->m_lower_bound;
            auto diff = lb - sum;
            auto coeff1 = cons->m_poly.m_var_factor_value[v1], coeff2 = cons->m_poly.m_var_factor_value[v2];
            if (cons->m_poly.m_vars.count(v1) != 0 && cons->m_poly.m_vars.count(v2) != 0)
            {
                cons->m_poly.store_factor_and_intercept_vars(v1);
                if (cons->m_poly.m_var_vars[v1].first.count(v2) != 0)
                { // Factor with each other
                    exist = false;
                    return;
                }
                // coeff1 * v1 + coeff2 * v2 + other == sum
                // coeff1 * (v1 + t_coeff1 * t) + coeff2 * (v2 + t_coeff2 * t) + other  >= lb
                // t_coeff * t + sum >= lb
                auto t_coeff = coeff1 * t_coeff1 + coeff2 * t_coeff2;
                if (t_coeff == 0)
                {
                    if (sum >= lb)
                    {
                        exist = true;
                        inf = true;
                    }
                    else
                    {
                        exist = false;
                    }
                }
                else if (t_coeff > 0)
                {
                    exist = true;
                    inf = true;
                }
                else
                {
                    exist = true;
                    inf = false;
                    val = std::floor(diff / t_coeff);
                }
            }
            else if (cons->m_poly.m_vars.count(v1) != 0 && cons->m_poly.m_vars.count(v2) == 0)
            {
                // coeff1 * v1 + other == sum
                // coeff1 * (v1 + t_coeff1 * t) + other >= lb
                // t_coeff * t + sum >= lb
                auto t_coeff = coeff1 * t_coeff1;
                if (t_coeff1 == 0)
                {
                    if (sum >= lb)
                    {
                        exist = true;
                        inf = true;
                    }
                    else
                    {
                        exist = false;
                    }
                }
                else if (t_coeff1 > 0)
                {
                    exist = true;
                    inf = true;
                }
                else
                {
                    exist = true;
                    inf = false;
                    val = std::floor(diff / t_coeff);
                }
            }
            else if (cons->m_poly.m_vars.count(v1) == 0 && cons->m_poly.m_vars.count(v2) != 0)
            {
                // coeff2 * v2 + other == sum
                // coeff2 * (v2 + t_coeff2 * t) + other >= lb
                // t_coeff * t + sum >= lb
                auto t_coeff = coeff2 * t_coeff2;
                if (t_coeff == 0)
                {
                    if (sum >= lb)
                    {
                        exist = true;
                        inf = true;
                    }
                    else
                    {
                        exist = false;
                    }
                }
                else if (t_coeff > 0)
                {
                    exist = true;
                    inf = true;
                }
                else
                {
                    exist = true;
                    inf = false;
                    val = std::floor(diff / t_coeff);
                }
            }
            else
            {
                exist = true;
                inf = true;
            }
        }

        void get_largest_vnd_feasible_int_cons_upper(constraint *cons, int v1, int v2, long double t_coeff1, long double t_coeff2, bool &exist, bool &inf, long double &val)
        {
            auto sum = cons->m_poly.get_sum();
            auto ub = cons->m_upper_bound;
            auto diff = ub - sum;
            auto coeff1 = cons->m_poly.m_var_factor_value[v1], coeff2 = cons->m_poly.m_var_factor_value[v2];
            if (cons->m_poly.m_vars.count(v1) != 0 && cons->m_poly.m_vars.count(v2) != 0)
            {
                cons->m_poly.store_factor_and_intercept_vars(v1);
                if (cons->m_poly.m_var_vars[v1].first.count(v2) != 0)
                { // Factor with each other
                    exist = false;
                    return;
                }
                // coeff1 * v1 + coeff2 * v2 + other == sum
                // coeff1 * (v1 + t_coeff1 * t) + coeff2 * (v2 + t_coeff2 * t) + other  <= ub
                // t_coeff * t + sum <= ub
                auto t_coeff = coeff1 * t_coeff1 + coeff2 * t_coeff2;
                if (t_coeff == 0)
                {
                    if (sum <= ub)
                    {
                        exist = true;
                        inf = true;
                    }
                    else
                    {
                        exist = false;
                    }
                }
                else if (t_coeff > 0)
                {
                    exist = true;
                    inf = false;
                    val = std::floor(diff / t_coeff);
                }
                else
                {
                    exist = true;
                    inf = true;
                }
            }
            else if (cons->m_poly.m_vars.count(v1) != 0 && cons->m_poly.m_vars.count(v2) == 0)
            {
                // coeff1 * v1 + other == sum
                // coeff1 * (v1 + t_coeff1 * t) + other <= ub
                // t_coeff * t + sum <= ub
                auto t_coeff = coeff1 * t_coeff1;
                if (t_coeff1 == 0)
                {
                    if (sum <= ub)
                    {
                        exist = true;
                        inf = true;
                    }
                    else
                    {
                        exist = false;
                    }
                }
                else if (t_coeff1 > 0)
                {
                    exist = true;
                    inf = false;
                    val = std::floor(diff / t_coeff);
                }
                else
                {
                    exist = true;
                    inf = true;
                }
            }
            else if (cons->m_poly.m_vars.count(v1) == 0 && cons->m_poly.m_vars.count(v2) != 0)
            {
                // coeff2 * v2 + other == sum
                // coeff2 * (v2 + t_coeff2 * t) + other <= ub
                // t_coeff * t + sum <= ub
                auto t_coeff = coeff2 * t_coeff2;
                if (t_coeff == 0)
                {
                    if (sum <= ub)
                    {
                        exist = true;
                        inf = true;
                    }
                    else
                    {
                        exist = false;
                    }
                }
                else if (t_coeff > 0)
                {
                    exist = true;
                    inf = false;
                    val = std::floor(diff / t_coeff);
                }
                else
                {
                    exist = true;
                    inf = true;
                }
            }
            else
            {
                exist = true;
                inf = true;
            }
        }

        /**
         * @brief Get the smallest int value in var's feasible set
         *
         * loop all constraints, get max of smallest value of each feasible set
         * lb(feasible_set) = max(lb(cons, v))
         *
         * 1. exist = e1/\e2/\...
         * 2. inf = inf1/\inf2/\...
         * 3. val = max(v1, v2, ...)
         */
        void get_smallest_feasible_int_value(int v, int &val, bool &inf, bool &exist)
        {
            exist = true;
            inf = true;
            val = INT32_MIN;
            // consider bound
            if (m_vars[v].m_lb)
            {
                inf = false;
                val = m_vars[v].m_lb;
            }
            // consider constraints
            for (auto idx : m_vars[v].m_constraints)
            {
                auto cons = m_constraints[idx];
                int curr_val;
                bool curr_inf, curr_exist;
                get_smallest_feasible_int_value_cons(v, cons, curr_val, curr_inf, curr_exist);
                if (!curr_exist)
                {
                    exist = false;
                    return;
                }
                inf = inf && curr_inf;
                if (!curr_inf)
                {
                    if (curr_val == val)
                    {
                        m_vars[v].m_crucial_constraint.push_back(idx);
                    }
                    else if (curr_val > val)
                    {
                        val = curr_val;
                        m_vars[v].m_crucial_constraint.clear();
                        m_vars[v].m_crucial_constraint.push_back(idx);
                    }
                }
            }
        }

        void get_smallest_feasible_int_value_cons(int v, constraint const *cons, int &val, bool &inf, bool &exist)
        {
            if (cons->m_lower && cons->m_upper)
            { // get larger one
                bool inf1, inf2, exist1, exist2;
                int val1, val2;
                get_smallest_feasible_int_value_cons_lower(v, cons, val1, inf1, exist1);
                get_smallest_feasible_int_value_cons_upper(v, cons, val2, inf2, exist2);
                if (!(exist1 && exist2))
                {
                    exist = false;
                }
                else
                {
                    exist = true;
                    if (inf1 && inf2)
                    {
                        inf = true;
                    }
                    else if (inf1 && !inf2)
                    {
                        val = val2;
                        inf = false;
                    }
                    else if (!inf1 && inf2)
                    {
                        val = val1;
                        inf = false;
                    }
                    else
                    {
                        inf = false;
                        val = std::max(val1, val2);
                    }
                }
            }
            else if (cons->m_lower && !cons->m_upper)
            {
                get_smallest_feasible_int_value_cons_lower(v, cons, val, inf, exist);
            }
            else if (!cons->m_lower && cons->m_upper)
            {
                get_smallest_feasible_int_value_cons_upper(v, cons, val, inf, exist);
            }
            else
            {
                UNREACHABLE();
            }
        }

        /**
         * a x + b >= c
         *
         * x >=< (c-b)/a  if a != 0
         */
        void get_smallest_feasible_int_value_cons_lower(int v, constraint const *cons, int &val, bool &inf, bool &exist)
        {
            long double sum = cons->m_poly.get_sum();
            auto coeff = cons->m_poly.m_var_factor_value.at(v);
            long double intercept = sum - coeff * m_assignment.value(v);
            if (coeff == 0)
            {
                if (sum >= cons->m_lower_bound)
                {
                    inf = true;
                    exist = true;
                }
                else
                {
                    exist = false;
                }
            }
            else
            {
                double goal = (cons->m_lower_bound - intercept) / coeff;
                if (coeff > 0)
                {
                    val = std::ceil(goal);
                    inf = false;
                    exist = true;
                }
                else
                {
                    inf = true;
                    exist = true;
                }
            }
        }

        /**
         * a x + b <= c
         *
         * x >=< (c-b)/a  if a != 0
         */
        void get_smallest_feasible_int_value_cons_upper(int v, constraint const *cons, int &val, bool &inf, bool &exist)
        {
            long double sum = cons->m_poly.get_sum();
            auto coeff = cons->m_poly.m_var_factor_value.at(v);
            long double intercept = sum - coeff * m_assignment.value(v);
            if (coeff == 0)
            {
                if (sum <= cons->m_upper_bound)
                {
                    exist = true;
                    inf = true;
                }
                else
                {
                    exist = false;
                }
            }
            else
            {
                double goal = (cons->m_upper_bound - intercept) / coeff;
                if (coeff > 0)
                {
                    inf = true;
                    exist = true;
                }
                else
                {
                    val = std::ceil(goal);
                    exist = true;
                    inf = false;
                }
            }
        }

        /**
         * @brief Get the largestest int value in var's feasible set
         *
         * loop all constraints, get min of smallest value of each feasible set
         * ub(feasible_set) = min(ub(cons, v))
         *
         * 1. exist = e1/\e2/\...
         * 2. inf = inf1/\inf2/\...
         * 3. val = min(v1, v2, ...)
         */
        void get_largest_feasible_int_value(int v, int &val, bool &inf, bool &exist)
        {
            exist = true;
            inf = true;
            val = INT32_MAX;
            // consider bound
            if (m_vars[v].m_lb)
            {
                inf = false;
                val = m_vars[v].m_lb;
            }
            // consider constraints
            for (auto idx : m_vars[v].m_constraints)
            {
                auto cons = m_constraints[idx];
                int curr_val;
                bool curr_inf, curr_exist;
                get_largest_feasible_int_value_cons(v, cons, curr_val, curr_inf, curr_exist);
                if (!curr_exist)
                {
                    exist = false;
                    return;
                }
                inf = inf && curr_inf;
                if (!curr_inf)
                {
                    if (curr_val == val)
                    {
                        m_vars[v].m_crucial_constraint.push_back(idx);
                    }
                    else if (curr_val < val)
                    {
                        val = curr_val;
                        m_vars[v].m_crucial_constraint.clear();
                        m_vars[v].m_crucial_constraint.push_back(idx);
                    }
                }
            }
        }

        void get_largest_feasible_int_value_cons(int v, constraint const *cons, int &val, bool &inf, bool &exist)
        {
            if (cons->m_lower && cons->m_upper)
            { // get larger one
                bool inf1, inf2, exist1, exist2;
                int val1, val2;
                get_largest_feasible_int_value_cons_lower(v, cons, val1, inf1, exist1);
                get_largest_feasible_int_value_cons_upper(v, cons, val2, inf2, exist2);
                if (!(exist1 && exist2))
                {
                    exist = false;
                }
                else
                {
                    exist = true;
                    if (inf1 && inf2)
                    {
                        inf = true;
                    }
                    else if (inf1 && !inf2)
                    {
                        val = val2;
                        inf = false;
                    }
                    else if (!inf1 && inf2)
                    {
                        val = val1;
                        inf = false;
                    }
                    else
                    {
                        inf = false;
                        val = std::min(val1, val2);
                    }
                }
            }
            else if (cons->m_lower && !cons->m_upper)
            {
                get_largest_feasible_int_value_cons_lower(v, cons, val, inf, exist);
            }
            else if (!cons->m_lower && cons->m_upper)
            {
                get_largest_feasible_int_value_cons_upper(v, cons, val, inf, exist);
            }
            else
            {
                UNREACHABLE();
            }
        }

        /**
         * a x + b >= c
         *
         * x >=< (c-b)/a  if a != 0
         */
        void get_largest_feasible_int_value_cons_lower(int v, constraint const *cons, int &val, bool &inf, bool &exist)
        {
            double sum = cons->m_poly.get_sum();
            auto coeff = cons->m_poly.m_var_factor_value.at(v);
            if (coeff == 0)
            {
                if (sum >= cons->m_lower_bound)
                {
                    inf = true;
                    exist = true;
                }
                else
                {
                    exist = false;
                }
            }
            else
            {
                double goal = (cons->m_lower_bound - sum) / coeff;
                if (coeff < 0)
                {
                    val = std::floor(goal);
                    inf = false;
                    exist = true;
                }
                else
                {
                    inf = true;
                    exist = true;
                }
            }
        }

        // Select vars meeting two prerequisites
        // 1. appear in objective
        // 2. is bounded by constraint not by its own bound
        int_table m_vars_in_obj_without_bound; // TODO: maintain

        void insert_and_select_paired_lift_moves(int &v1, int &v2, long double &shift1, long double &shift2, long double &score)
        {
            // Step 1. Select var in obj without bound
            m_vars_in_obj_without_bound.clear();
            for (int v : m_vars_in_obj)
            {
                auto val = m_assignment.value(v);
                if (m_vars[v].m_coeff > 0)
                {
                    if (!(m_vars[v].m_lower && val == m_vars[v].m_lb))
                    {
                        m_vars_in_obj_without_bound.insert(v);
                    }
                }
                else if (m_vars[v].m_coeff < 0)
                {
                    if (!(m_vars[v].m_upper && val == m_vars[v].m_ub))
                    {
                        m_vars_in_obj_without_bound.insert(v);
                    }
                }
            }
            // Step 2. Insert var pairs
            int_pair_vector m_operated_vars;
            std::vector<direction> m_operated_dirs;
            // Record whether other_v has been inserted for each obj_v
            int_table curr_vars;
            std::vector<int> changed_vars;
            for (int obj_v : m_vars_in_obj_without_bound)
            {
                // ASSERT(!m_vars[obj_v].m_crucial_constraint.empty());
                curr_vars.clear();
                // if(m_assignment.value(obj_v) ==0) continue;
                for (auto idx : m_vars[obj_v].m_constraints)
                {
                    auto cons = m_constraints[idx];
                    cons->m_poly.store_factor_and_intercept_vars(obj_v);
                    ASSERT(cons->m_poly.m_var_factor_value.count(obj_v) != 0);
                    auto coeff1 = cons->m_poly.m_var_factor_value.at(obj_v);
                    // Ignore coeff1 == 0 case
                    // ASSERT(coeff1 != 0);
                    if (coeff1 == 0)
                        continue;
                    // Loop constraint's poly's intercept at var v
                    ASSERT(cons->m_poly.m_var_vars.count(obj_v) != 0);
                    // for(int other_v: cons->m_poly.m_var_vars.at(obj_v).first)
                    // {
                    //     if(m_assignment.value(other_v) !=0) changed_vars.push_back(other_v);
                    // }
                    for (int other_v : cons->m_poly.m_var_vars.at(obj_v).second)
                    {
                        // Vars in obj, inserted already or also appear in factor will not be considered
                        if (m_vars_in_obj.count(other_v) != 0 || curr_vars.count(other_v) != 0 || cons->m_poly.m_var_vars.at(obj_v).first.count(other_v) != 0)
                        {
                            continue;
                        }
                        ASSERT(cons->m_poly.m_var_factor_value.count(other_v) != 0);
                        auto coeff2 = cons->m_poly.m_var_factor_value.at(other_v);
                        // Ignore coeff2 == 0 case
                        if (coeff2 == 0)
                        {
                            continue;
                        }
                        auto coeff_ratio = -1.0 * coeff2 / coeff1;
                        // coeff1 * obj_v + coeff2 * other_v >= lb
                        if (cons->m_lower)
                        {
                            if (m_vars[obj_v].m_coeff > 0)
                            { // we make change to let obj_v get smaller
                                if (coeff1 > 0)
                                { // obj_v >= (lb - coeff2 * other_v) / coeff1
                                    if (coeff_ratio > 0)
                                    { // we make other_v smaller
                                        curr_vars.insert(other_v);
                                        m_operated_vars.push_back(std::make_pair(other_v, obj_v));
                                        m_operated_dirs.push_back(down);
                                    }
                                    else
                                    { // we make other_v larger
                                        curr_vars.insert(other_v);
                                        m_operated_vars.push_back(std::make_pair(other_v, obj_v));
                                        m_operated_dirs.push_back(up);
                                    }
                                }
                                else
                                { // obj_v <= goal, do not consider
                                }
                            }
                            else
                            { // we make change to let obj_v get larger
                                if (coeff1 > 0)
                                { // obj_v >= goal, do not consider
                                }
                                else
                                { // obj_v <= (lb - coeff2 * other_v) / coeff_1
                                    if (coeff_ratio > 0)
                                    { // we make other_v larger
                                        curr_vars.insert(other_v);
                                        m_operated_vars.push_back(std::make_pair(other_v, obj_v));
                                        m_operated_dirs.push_back(up);
                                    }
                                    else
                                    { // we make other_v smaller
                                        curr_vars.insert(other_v);
                                        m_operated_vars.push_back(std::make_pair(other_v, obj_v));
                                        m_operated_dirs.push_back(down);
                                    }
                                }
                            }
                        }
                        // coeff1 * obj_v + coeff2 * other_v <= ub
                        if (cons->m_upper)
                        {
                            if (m_vars[obj_v].m_coeff > 0)
                            { // we make change to let obj_v smaller
                                if (coeff1 > 0)
                                { // obj <= goal, do not consider
                                }
                                else
                                { // obj >= (ub - coeff2 * other_v) / coeff1
                                    if (coeff_ratio > 0)
                                    { // we make other_v smaller
                                        curr_vars.insert(other_v);
                                        m_operated_vars.push_back(std::make_pair(other_v, obj_v));
                                        m_operated_dirs.push_back(down);
                                    }
                                    else
                                    { // we make other_v larger
                                        curr_vars.insert(other_v);
                                        m_operated_vars.push_back(std::make_pair(other_v, obj_v));
                                        m_operated_dirs.push_back(up);
                                    }
                                }
                            }
                            else
                            { // we make change to let obj_v larger
                                if (coeff1 > 0)
                                { // obj <= (ub - coeff2 * other_v) / coeff1
                                    if (coeff_ratio > 0)
                                    { // we make other_v larger
                                        curr_vars.insert(other_v);
                                        m_operated_vars.push_back(std::make_pair(other_v, obj_v));
                                        m_operated_dirs.push_back(up);
                                    }
                                    else
                                    { // we make other_v smaller
                                        curr_vars.insert(other_v);
                                        m_operated_vars.push_back(std::make_pair(other_v, obj_v));
                                        m_operated_dirs.push_back(down);
                                    }
                                }
                                else
                                { // obj >= goal, do not consider
                                }
                            }
                        }
                    }
                }
            }
            // Step 3. Calculate shift
            bool BMS;
            int cnt;
            // std::cout<<" pair size: "<< m_operated_vars.size()<<std::endl;
            if (m_operated_vars.size() >= DOUBLE_BMS_VAL)
            {
                BMS = true;
                cnt = DOUBLE_BMS_VAL;
            }
            else
            {
                BMS = false;
                cnt = m_operated_vars.size();
            }
            int size = m_operated_vars.size();
            int obj_v, other_v;
            long double obj_shift, other_shift, curr_score;
            direction curr_dir;
            score = INT32_MIN;
            m_var_large_threshold.clear();
            m_var_small_threshold.clear();
            for (int i = 0; i < cnt; i++)
            {
                if (BMS)
                {
                    int rand_index = rand() % (size - i);
                    other_v = m_operated_vars[rand_index].first;
                    obj_v = m_operated_vars[rand_index].second;
                    curr_dir = m_operated_dirs[rand_index];
                    m_operated_vars[rand_index] = m_operated_vars[size - i - 1];
                    m_operated_dirs[rand_index] = m_operated_dirs[size - i - 1];
                }
                else
                {
                    other_v = m_operated_vars[i].first;
                    obj_v = m_operated_vars[i].second;
                    curr_dir = m_operated_dirs[i];
                }
                calculate_pair_lift_shift(other_v, obj_v, other_shift, obj_shift, curr_score, curr_dir);
                if (curr_score > score)
                {
                    v1 = other_v;
                    v2 = obj_v;
                    shift1 = other_shift;
                    shift2 = obj_shift;
                    score = curr_score;
                }
            }
        }

        std::unordered_map<int, int> m_var_large_threshold, m_var_small_threshold;

        /*
        v1: first var, or other var
        v2: second var, or lifted var
        */
        void calculate_pair_lift_shift(int other_v, int obj_v, long double &other_shift, long double &obj_shift, long double &score, direction dir)
        {
            score = INT32_MIN;
            // Step I. Calculate shift of other_v (threshold value in feasible set)
            if (dir == up)
            {
                if (m_var_large_threshold.count(other_v) == 0)
                {
                    bool inf, exist;
                    int val;
                    get_largest_feasible_int_value(other_v, val, inf, exist);
                    if (exist && !inf)
                    {
                        other_shift = val - m_assignment.value(other_v);
                        m_var_large_threshold[other_v] = other_shift;
                    }
                }
                else
                {
                    other_shift = m_var_large_threshold[other_v];
                }
            }
            else
            {
                if (m_var_small_threshold.count(other_v) == 0)
                {
                    bool inf, exist;
                    int val;
                    get_smallest_feasible_int_value(other_v, val, inf, exist);
                    if (exist && !inf)
                    {
                        other_shift = val - m_assignment.value(other_v);
                        m_var_small_threshold[other_v] = other_shift;
                    }
                    else
                    {
                        return;
                    }
                }
                else
                {
                    other_shift = m_var_small_threshold[other_v];
                }
            }
            // Step II. Given a shift for other_v, calculate shift for obj_v
            m_assignment.shift(other_v, other_shift);
            long_double_vector m_sum_cached, m_factor_cached;
            for (auto idx : m_vars[other_v].m_constraints)
            {
                auto cons = m_constraints[idx];
                if (cons->m_poly.m_vars.count(obj_v) == 0)
                {
                    continue;
                }
                auto pre_sum = cons->m_poly.get_sum();
                m_sum_cached.push_back(pre_sum);
                auto post_sum = pre_sum + cons->m_poly.m_var_factor_value[other_v] * other_shift;
                cons->m_poly.set_sum(post_sum);
                m_factor_cached.push_back(cons->m_poly.m_var_factor_value[obj_v]);
                adjust_factor_value(obj_v, cons);
            }
            if (m_vars[obj_v].m_coeff > 0)
            {
                bool exist, inf;
                int val;
                get_smallest_feasible_int_value(obj_v, val, inf, exist);
                if (exist && !inf)
                {
                    obj_shift = val - m_assignment.value(obj_v);
                    score = -1.0 * m_vars[obj_v].m_coeff * obj_shift;
                }
                else
                {
                    score = INT32_MIN;
                }
            }
            else
            {
                bool exist, inf;
                int val;
                get_largest_feasible_int_value(obj_v, val, inf, exist);
                if (exist && !inf)
                {
                    obj_shift = val - m_assignment.value(obj_v);
                    score = -1.0 * m_vars[obj_v].m_coeff * obj_shift;
                }
                else
                {
                    score = INT32_MIN;
                }
            }
            // Step III. Rstore original status for constraints
            m_assignment.shift(other_v, -1.0 * other_shift);
            int i = 0;
            for (auto idx : m_vars[other_v].m_constraints)
            {
                auto cons = m_constraints[idx];
                if (cons->m_poly.m_vars.count(obj_v) == 0)
                {
                    continue;
                }
                cons->m_poly.set_sum(m_sum_cached[i]);
                cons->m_poly.m_var_factor_value[obj_v] = m_factor_cached[i];
                i++;
            }
        }

        void execute_pair_move(int v1, int v2, long double shift1, long double shift2)
        {
            execute_critical_move(v1, shift1);
            execute_critical_move(v2, shift2);
        }

        /**
         * a x + b <= c
         *
         * x >=< (c-b)/a  if a != 0
         */
        void get_largest_feasible_int_value_cons_upper(int v, constraint const *cons, int &val, bool &inf, bool &exist)
        {
            double sum = cons->m_poly.get_sum();
            auto coeff = cons->m_poly.m_var_factor_value.at(v);
            if (coeff == 0)
            {
                if (sum <= cons->m_upper_bound)
                {
                    exist = true;
                    inf = true;
                }
                else
                {
                    exist = false;
                }
            }
            else
            {
                double goal = (cons->m_upper_bound - sum) / coeff;
                if (coeff < 0)
                {
                    inf = true;
                    exist = true;
                }
                else
                {
                    val = std::floor(goal);
                    exist = true;
                    inf = false;
                }
            }
        }

        void execute_vnd_move(int v1, int v2, long double shift1, long double shift2)
        {
            execute_critical_move(v1, shift1);
            update_best_information();
            execute_critical_move(v2, shift2);
        }

        static inline bool is_int(long double x)
        {
            return std::floor(x) == x;
        }

        void solve()
        {
            enable_switch_choices();
            TRACE(display_objective(tout);
                  display_all_constraints_with_weight(tout););
            m_restart = 0;
            m_lift_step = 0;
            m_lift_pair_step = 0;
            m_vnd_lift_step = 0;
            std::cout << "#vars: " << m_vars.size() << std::endl;
            std::cout << "#demand: " << m_demand_cons.size() << std::endl;
            std::cout << "#cons: " << m_constraints.size() << std::endl;
            initialize();
            TRACE(
                tout << "after initialize\n";
                display_assignment(tout, m_assignment););
            for (m_step = 0; m_step <= max_step; m_step++)
            {
                std::cout << "step: " << m_step << ", ";
                std::cout << "#unsat: " << m_unsat_constraints.size() << std::endl;
                TRACE(
                    tout << "step: " << m_step << std::endl;
                    display_unsat_constraints(tout);
                    display_all_constraints_with_weight(tout););
                // check_sat_constraints();
                // Succeed
                if ((enable_cutoff && time_elapsed() >= m_cutoff && m_step % 100 == 0) || (enable_stepoff && m_step >= m_stepoff))
                {
                    std::cout << "#step: " << m_step << std::endl;
                    std::cout << "#time: " << time_elapsed() << std::endl;
                    TRACE(display_find_sat(tout););
                    check_result();
                    std::cout << "------------ pass ------------";
                    display_results(std::cout);
                    return;
                }
                // Main Search
                // sat status
                if (is_sat_status)
                {
                    if (!enable_all_lift_move_pool)
                    {
                        int var_index = INT32_MAX;
                        long double score, var_shift = INT32_MAX;
                        if (enable_lift_move)
                        {
                            // if (enable_vnd_lift_move) {
                            //     int v1, v2;
                            //     long double shift1, shift2;
                            //     insert_and_select_vnd_lift_moves(v1, v2, shift1, shift2, score);
                            //     if (score > 0) {
                            //         m_vnd_lift_step++;
                            //         execute_vnd_move(v1, v2, shift1, shift2);
                            //         update_best_information();
                            //         continue;
                            //     }
                            // }
                            insert_and_select_lift_moves(var_index, var_shift, score);
                            if (score > 0)
                            {
                                std::cout << "lift " << std::endl;
                                m_lift_step++;
                                execute_critical_move(var_index, var_shift);
                                update_best_information();
                                continue;
                            }
                            // if(enable_paired_lift_move) {
                            //     int v1, v2;
                            //     long double shift1, shift2;
                            //     insert_and_select_paired_lift_moves(v1, v2, shift1, shift2, score);
                            //     if(score > 0) {
                            //         m_lift_pair_step++;
                            //         execute_pair_move(v1, v2, shift1, shift2);
                            //         update_best_information();
                            //         continue;
                            //     }
                            // }
                        }
                        insert_operation_from_unbound_constraints_objective();
                        TRACE(display_operation_table(tout, m_operation_table););
                        select_best_operation(var_index, var_shift, score);
                        TRACE(display_best_choice(tout, var_index, var_shift, score););
                        // Greedy
                        if (score > 0)
                        {
                            execute_critical_move(var_index, var_shift);
                        }
                        else
                        { // Random
                            update_weight();
                            TRACE(tout << "after weighting\n";
                                  display_all_constraints_with_weight(tout);
                                  tout << "object weight: " << m_object_weight << std::endl;);
                            random_walk_sat();
                        }
                    }
                    else
                    { // Put life move, paired lift move and vnd lift move in a pool

                        // new
                        //  int var_index = INT32_MAX;
                        //  long double score, var_shift = INT32_MAX;
                        //  insert_and_select_lift_moves(var_index, var_shift, score);
                        //  if(score > 0) {
                        //      // std::cout <<"lift "<<std::endl;
                        //      m_lift_step++;
                        //      execute_critical_move(var_index, var_shift);
                        //      update_best_information();
                        //      continue;
                        //  }
                        // new
                        int single_var = INT32_MAX;
                        long double single_score, single_shift = INT32_MAX;
                        insert_and_select_lift_moves(single_var, single_shift, single_score);
                        int vnd_v1, vnd_v2;
                        long double vnd_shift1, vnd_shift2, vnd_score;
                        insert_and_select_vnd_lift_moves(vnd_v1, vnd_v2, vnd_shift1, vnd_shift2, vnd_score);
                        int pair_v1, pair_v2;
                        long double pair_shift1, pair_shift2, pair_score;
                        insert_and_select_paired_lift_moves(pair_v1, pair_v2, pair_shift1, pair_shift2, pair_score);
                        vnd_score = INT32_MIN;
                        // pair_score = INT32_MIN;
                        // single_score = INT32_MIN;
                        if (single_score <= 0 && vnd_score <= 0 && pair_score <= 0)
                        { // no positive score exists
                            int var_index;
                            long double score, var_shift;
                            insert_operation_from_unbound_constraints_objective();
                            TRACE(display_operation_table(tout, m_operation_table););
                            select_best_operation(var_index, var_shift, score);
                            TRACE(display_best_choice(tout, var_index, var_shift, score););
                            // Greedy
                            if (score > 0)
                            {
                                execute_critical_move(var_index, var_shift);
                            }
                            else
                            { // Random
                                update_weight();
                                TRACE(tout << "after weighting\n";
                                      display_all_constraints_with_weight(tout);
                                      tout << "object weight: " << m_object_weight << std::endl;);
                                random_walk_sat();
                            }
                        }
                        else
                        {
                            long double best_score = std::max(single_score, std::max(vnd_score, pair_score));
                            if (best_score == single_score)
                            {
                                std::cout << "singel " << std::endl;
                                std::cout << single_score << pair_score << std::endl;
                                execute_critical_move(single_var, single_shift);
                            }
                            else if (best_score == vnd_score)
                            {
                                std::cout << "vnd " << std::endl;
                                execute_vnd_move(vnd_v1, vnd_v2, vnd_shift1, vnd_shift2);
                            }
                            else if (best_score == pair_score)
                            {
                                std::cout << "pair " << std::endl;
                                execute_pair_move(pair_v1, pair_v2, pair_shift1, pair_shift2);
                            }
                            else
                            {
                                UNREACHABLE();
                            }
                        }
                    }
                }
                else
                { // unsat status
                    if (enable_vnd_search)
                    { // enable vnd version
                      // insert_operation_from_unsat_constraints(); // single move
                      // TRACE(display_operation_table(tout, m_operation_table););
                      // int var_index = INT32_MAX;
                      // long double score, var_shift = INT32_MAX;
                      // select_best_operation(var_index, var_shift, score); // select best single move
                      // TRACE(display_best_choice(tout, var_index, var_shift, score););
                      // if(score > 0) { // greedy single move
                      //     execute_critical_move(var_index, var_shift); // execute single move
                      // } else {
                      //     insert_vnd_move_from_unsat_level_one(); // only from intercept
                      //     int v1, v2; long double shift1, shift2;
                      //     select_best_vnd_move(v1, v2, shift1, shift2, score); // select best vnd move
                      //     if(score > 0) { // greedy vnd level one
                      //         m_vnd_step++;
                      //         execute_vnd_move(v1, v2, shift1, shift2);
                      //     } else {
                      //         if(enable_vnd_two_level) {
                      //             m_vnd_step_two_level++;
                      //             insert_vnd_move_from_unsat_level_two(); // only from factor
                      //             select_best_vnd_move(v1, v2, shift1, shift2, score);
                      //             if (score > 0) {
                      //                 execute_vnd_move(v1, v2, shift1, shift2);
                      //                 update_best_information();
                      //                 continue;
                      //             }
                      //         }
                      //         // Random walk
                      //         update_weight();
                      //         TRACE(tout << "after weighting\n";
                      //                 display_all_constraints_with_weight(tout);
                      //                 tout << "object weight: " << m_object_weight << std::endl;);
                      //         random_walk_unsat();
                      //     }
                      // }
                    }
                    else
                    { // original version
                        insert_operation_from_unsat_constraints();
                        TRACE(display_operation_table(tout, m_operation_table););
                        int var_index = INT32_MAX;
                        long double score, var_shift = INT32_MAX;
                        select_best_operation(var_index, var_shift, score, true);
                        TRACE(display_best_choice(tout, var_index, var_shift, score););
                        // Greedy
                        if (score > 0)
                        {
                            execute_critical_move(var_index, var_shift);
                        }
                        else
                        { // Random
                            update_weight();
                            TRACE(tout << "after weighting\n";
                                  display_all_constraints_with_weight(tout);
                                  tout << "object weight: " << m_object_weight << std::endl;);
                            random_walk_unsat();
                        }
                    }
                }
                update_best_information();
            }
        }

        /*
            + 3 x1
            + [3 x2 * x3]
            - [4 x1 * x3]
        */
        std::ostream &display_monomial(std::ostream &out, monomial const &mono) const
        {
            if (mono.m_coeff == 0 || mono.size() == 0)
            {
                return out;
            }
            if (mono.size() == 1)
            {
                int var_idx = *(mono.m_vars.begin());
                if (mono.m_coeff > 0)
                {
                    out << " +";
                }
                out << " " << mono.m_coeff << " " << m_vars[var_idx].m_name;
            }
            else
            {
                if (mono.m_coeff > 0)
                {
                    out << " +";
                }
                out << " " << mono.m_coeff << " ";
                int index = 0;
                for (int var_idx : mono.m_vars)
                {
                    out << m_vars[var_idx].m_name;
                    if (index != mono.m_vars.size() - 1)
                    {
                        out << " * ";
                    }
                    index++;
                }
            }
            return out;
        }

        /*
            3 x1
            [ -3 x2 * x3 ]
        */
        std::ostream &display_first_monomial(std::ostream &out, monomial const &mono) const
        {
            if (mono.m_coeff == 0 || mono.size() == 0)
            {
                return out;
            }
            if (mono.size() == 1)
            {
                int var_idx = *(mono.m_vars.begin());
                out << " " << mono.m_coeff << " " << m_vars[var_idx].m_name;
            }
            else
            {
                int index = 0;
                out << mono.m_coeff << " ";
                for (int var_idx : mono.m_vars)
                {
                    out << m_vars[var_idx].m_name;
                    if (index != mono.m_vars.size() - 1)
                    {
                        out << " * ";
                    }
                    index++;
                }
            }
            return out;
        }

        std::ostream &display_polynomial(std::ostream &out, polynomial const &poly) const
        {
            int index = 0;
            // display linear part
            for (auto idx : poly.m_linear_terms)
            {
                if (index == 0)
                {
                    display_first_monomial(out, poly.m_monomials[idx]);
                }
                else
                {
                    display_monomial(out, poly.m_monomials[idx]);
                }
                index++;
            }
            // display nonlinear part
            if (poly.m_nonlinear_terms.empty())
            {
                return out;
            }
            if (!poly.m_linear_terms.empty())
            {
                out << " +";
            }
            out << " [ ";
            index = 0;
            for (auto idx : poly.m_nonlinear_terms)
            {
                if (index == 0)
                {
                    display_first_monomial(out, poly.m_monomials[idx]);
                }
                else
                {
                    display_monomial(out, poly.m_monomials[idx]);
                }
                index++;
            }
            // hexiang
            out << " ]";
            return out;
        }

        std::ostream &display_constraint(std::ostream &out, constraint const &cons) const
        {
            if (cons.m_lower)
            {
                display_constraint_lower(out, cons);
            }
            if (cons.m_upper)
            {
                display_constraint_upper(out, cons);
            }
            return out;
        }

        std::ostream &display_constraint_with_weight(std::ostream &out, constraint const &cons) const
        {
            if (cons.m_lower)
            {
                display_constraint_lower_with_weight(out, cons);
            }
            if (cons.m_upper)
            {
                display_constraint_upper_with_weight(out, cons);
            }
            return out;
        }

        std::ostream &display_constraint_lower(std::ostream &out, constraint const &cons) const
        {
            out << cons.m_name + "_lower: ";
            display_polynomial(out, cons.m_poly);
            out << " >= " << cons.m_lower_bound << std::endl;
            return out;
        }

        std::ostream &display_constraint_upper(std::ostream &out, constraint const &cons) const
        {
            out << cons.m_name + "_upper: ";
            display_polynomial(out, cons.m_poly);
            out << " <= " << cons.m_upper_bound << std::endl;
            return out;
        }

        std::ostream &display_constraint_lower_with_weight(std::ostream &out, constraint const &cons) const
        {
            out << cons.m_name + "_lower (" << cons.m_lower_weight << "): ";
            display_polynomial(out, cons.m_poly);
            out << " >= " << cons.m_lower_bound << std::endl;
            return out;
        }

        std::ostream &display_constraint_upper_with_weight(std::ostream &out, constraint const &cons) const
        {
            out << cons.m_name + "_upper (" << cons.m_upper_weight << "): ";
            display_polynomial(out, cons.m_poly);
            out << " <= " << cons.m_upper_bound << std::endl;
            return out;
        }

        std::ostream &display_var_bound(std::ostream &out, var const &v) const
        {
            if (!v.m_lower && !v.m_upper)
            {
                return out;
            }
            else if (v.m_lower && !v.m_upper)
            {
                out << " " << v.m_lb << " <= " << v.m_name << std::endl;
            }
            else if (!v.m_lower && v.m_upper)
            {
                out << " " << v.m_name << " <= " << v.m_ub << std::endl;
            }
            else
            {
                out << " " << v.m_lb << " <= " << v.m_name << " <= " << v.m_ub << std::endl;
            }
            return out;
        }

        std::ostream &display_objective(std::ostream &out) const
        {
            monomial_vector m_monomials;
            int_table vars;
            for (auto v : m_vars)
            {
                if (v.m_coeff != 0)
                {
                    vars.clear();
                    vars.insert(v.m_index);
                    m_monomials.push_back(monomial(v.m_coeff, vars));
                }
            }
            display_polynomial(out, polynomial(m_monomials));
            out << std::endl;
            return out;
        }

        std::ostream &display_results(std::ostream &out) const
        {
            out << "-------------- RESULT --------------\n";
            out << "#step: " << m_step << std::endl;
            out << "#time: " << time_elapsed() << std::endl;
            out << "#lift step: " << m_lift_step << std::endl;
            out << "#lift pair step: " << m_lift_pair_step << std::endl;
            out << "#vnd lift step: " << m_vnd_lift_step << std::endl;
            out << "#vnd step: " << m_vnd_step << std::endl;
            out << "#vnd level two step: " << m_vnd_step_two_level << std::endl;
            display_find_sat(out);
            return out;
        }

        // void write_lp_file(std::string str)
        // {
        //     for (int i = 0; i < 20; i++)
        //     {
        //         std::string pre = "/pub/netdisk1/lijy/alimama/lp_learning/";
        //         std::string filename = pre + "test" + "_" + std::to_string(i) + ".lp";
        //         std::cout << filename << std::endl;
        //         std::ofstream outFile(filename);
        //         outFile << "Minimize\n";
        //         outFile << " obj: ";
        //         display_objective(outFile) << std::endl;
        //         outFile << "Subject To\n";
        //         initialize();
        //         int rand_index = rand() % m_vars_in_obj.size();
        //         var var_cur = m_vars[rand_index];
        //         std::unordered_set<int> used_vars;
        //         for (int used_var : m_vars_in_obj)
        //         {
        //             used_vars.insert(used_var);
        //         }
        //         for (auto ele : var_cur.m_constraints)
        //         {
        //             auto cons = m_constraints[ele];
        //             for (int used_var : cons->m_poly.m_vars)
        //             {
        //                 used_vars.insert(used_var);
        //             }
        //             display_constraint(outFile, *cons);
        //         }
        //         outFile << "\n";
        //         outFile << "Bounds\n";
        //         for (auto v : used_vars)
        //         {
        //             var cur = m_vars[v];
        //             display_var_bound(outFile, cur);
        //         }
        //         outFile << "\n";
        //         outFile << "Generals\n";
        //         for (auto v : used_vars)
        //         {
        //             var cur = m_vars[v];
        //             outFile << " " << cur.m_name;
        //         }
        //         outFile << "\n";
        //         outFile.close();
        //     }
        // }

        void write_lp_file(std::string str) {
            std::ofstream outFile(str);
            outFile << "Minimize\n";
            outFile << " obj: ";
            display_objective(outFile) << std::endl;
            outFile << "Subject To\n";
            for(auto ele: m_constraints) {
                display_constraint(outFile, *ele);
            }
            outFile << "\n";
            outFile << "Bounds\n";
            for(auto v: m_vars) {
                display_var_bound(outFile, v);
            }
            outFile << "\n";
            outFile << "Generals\n";
            for(auto v: m_vars) {
                outFile << " " << v.m_name;
            }
            outFile << "\n";            
            outFile << "End\n";
            outFile.close();
        }
        
        static bool is_finite_number(long double b)
        {
            return b >= -__DBL_MAX__ && b <= __DBL_MAX__;
        }

        static bool is_number(long double b)
        {
            return b == b;
        }
    };

    opt_solver::opt_solver()
    {
        m_imp = new imp();
    }

    void opt_solver::set_cutoff(long double _cutoff)
    {
        m_imp->set_cutoff(_cutoff);
    }

    void opt_solver::set_stepoff(int _step)
    {
        m_imp->set_stepoff(_step);
    }

    void opt_solver::enable_greedy_initialize()
    {
        m_imp->enable_greedy_initialize();
    }

    int opt_solver::register_var(std::string str, long double _coeff, bool lower, bool upper, long double lb, long double ub)
    {
        return m_imp->register_var(str, _coeff, lower, upper, lb, ub);
    }

    void opt_solver::register_constraint(polynomial const &poly, constraint_kind _kind, std::string _name, bool lower, bool upper, long double lb, long double ub)
    {
        m_imp->register_constraint(poly, _kind, _name, lower, upper, lb, ub);
    }

    void opt_solver::solve()
    {
        m_imp->solve();
    }

    void opt_solver::write_lp_file(std::string str)
    {
        m_imp->write_lp_file(str);
    }

    std::ostream &opt_solver::display_results(std::ostream &out) const
    {
        return m_imp->display_results(out);
    }
};