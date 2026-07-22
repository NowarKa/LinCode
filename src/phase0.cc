#include "phase0.hh"

#include <memory>
#include <scip/scip.h>
#include <scip/scip_var.h>
#include <scip/scipdefplugins.h>
#include <vector>

using namespace std;

auto check_feasibility(const Equations &system, int r,
                       shared_ptr<long long int> time_scip) -> bool
{

  chrono::steady_clock::time_point begin = chrono::steady_clock::now();

  // Initialization
  SCIP *scip = nullptr;
  SCIPcreate(&scip);
  SCIPincludeDefaultPlugins(scip);
  SCIPcreateProbBasic(scip, "phase0");

  // Variables
  vector<SCIP_VAR *> vars(system.A[0].size());

  for (size_t j = 0; j < vars.size(); j++)
  {
    SCIPcreateVarBasic(scip, &vars[j], ("x" + to_string(j)).c_str(), 0.0,
                       system.upper_bounds[j], 0.0, SCIP_VARTYPE_INTEGER);

    SCIPaddVar(scip, vars[j]);
  }

  // Adding bound (for r >= 2)
  auto np_kp1 = system.A[0].size() / 2;
  vector<SCIP_VAR *> u(system.A[0].size() / 2);

  for (size_t p = 0; p < np_kp1; p++)
  {
    SCIPcreateVarBasic(scip, &u[p], ("u" + to_string(p)).c_str(), 0.0, 1.0, 0.0,
                       SCIP_VARTYPE_BINARY);

    SCIPaddVar(scip, u[p]);
  }

  // xP <= ΛP * uP
  for (size_t p = 0; p < np_kp1; ++p)
  {
    SCIP_CONS *cons;

    SCIPcreateConsBasicLinear(scip, &cons, ("upper_" + to_string(p)).c_str(), 0,
                              nullptr, nullptr, -SCIPinfinity(scip), 0.0);

    SCIPaddCoefLinear(scip, cons, vars[p], 1.0);

    SCIPaddCoefLinear(scip, cons, u[p], -system.upper_bounds[p]);

    SCIPaddCons(scip, cons);
    SCIPreleaseCons(scip, &cons);
  }

  // xp >= r * up
  for (size_t p = 0; p < np_kp1; ++p)
  {
    SCIP_CONS *cons;

    SCIPcreateConsBasicLinear(scip, &cons, ("lower_" + to_string(p)).c_str(), 0,
                              nullptr, nullptr, 0.0, SCIPinfinity(scip));

    SCIPaddCoefLinear(scip, cons, vars[p], 1.0);

    SCIPaddCoefLinear(scip, cons, u[p], -static_cast<double>(r));

    SCIPaddCons(scip, cons);
    SCIPreleaseCons(scip, &cons);
  }

  // Ax = d
  for (size_t i = 0; i < system.A.size(); ++i)
  {
    SCIP_CONS *cons;

    SCIPcreateConsBasicLinear(scip, &cons, ("c" + to_string(i)).c_str(), 0,
                              nullptr, nullptr, system.d[i], system.d[i]);
    for (size_t j = 0; j < vars.size(); ++j)
    {
      if (system.A[i][j] != 0)
      {
        SCIPaddCoefLinear(scip, cons, vars[j], system.A[i][j]);
      }
    }
    SCIPaddCons(scip, cons);
    SCIPreleaseCons(scip, &cons);
  }

  // Check feasible only: one solution
  SCIPsetIntParam(scip, "limits/solutions", 1);

  // Resolution
  SCIPsolve(scip);
  SCIP_STATUS status = SCIPgetStatus(scip);
  bool feasible =
      status == SCIP_STATUS_OPTIMAL || status != SCIP_STATUS_INFEASIBLE;

  // Free memory
  for (auto &v : vars)
    SCIPreleaseVar(scip, &v);

  for (auto &v : u)
    SCIPreleaseVar(scip, &v);

  SCIPfree(&scip);

  chrono::steady_clock::time_point end = chrono::steady_clock::now();

  *time_scip +=
      chrono::duration_cast<chrono::milliseconds>(end - begin).count();

  return feasible;
}
