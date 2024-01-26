#include "slothApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
slothApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

slothApp::slothApp(InputParameters parameters) : MooseApp(parameters)
{
  slothApp::registerAll(_factory, _action_factory, _syntax);
}

slothApp::~slothApp() {}

void 
slothApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<slothApp>(f, af, s);
  Registry::registerObjectsTo(f, {"slothApp"});
  Registry::registerActionsTo(af, {"slothApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
slothApp::registerApps()
{
  registerApp(slothApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
slothApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  slothApp::registerAll(f, af, s);
}
extern "C" void
slothApp__registerApps()
{
  slothApp::registerApps();
}
