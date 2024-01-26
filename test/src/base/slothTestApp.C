//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "slothTestApp.h"
#include "slothApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
slothTestApp::validParams()
{
  InputParameters params = slothApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

slothTestApp::slothTestApp(InputParameters parameters) : MooseApp(parameters)
{
  slothTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

slothTestApp::~slothTestApp() {}

void
slothTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  slothApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"slothTestApp"});
    Registry::registerActionsTo(af, {"slothTestApp"});
  }
}

void
slothTestApp::registerApps()
{
  registerApp(slothApp);
  registerApp(slothTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
slothTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  slothTestApp::registerAll(f, af, s);
}
extern "C" void
slothTestApp__registerApps()
{
  slothTestApp::registerApps();
}
