sub updateGlobalVariable(oProject, varName, varValue)
oProject.ChangeProperty Array("NAME:AllTabs", Array("NAME:ProjectVariableTab", Array("NAME:PropServers",  _
  "ProjectVariables"), Array("NAME:ChangedProps", Array("NAME:$"+varName, "Value:=",  _
  varValue))))
end sub

Dim oAnsoftApp
Dim oDesktop
Dim oProject
Dim oDesign
Dim oEditor
Dim oModule
Set oAnsoftApp = CreateObject("Ansoft.ElectronicsDesktop")
Set oDesktop = oAnsoftApp.GetAppDesktop()
oDesktop.RestoreWindow
Set oProject = oDesktop.GetActiveProject()


call updateGlobalVariable(oProject, "x_pts", "[0,28.2900600918,29.089779105,29.6421998266,42.6336904312,43.2,48,46.3644396619,27.335700884,28.3,29.8470569183,19.3,24.8849430626,27.5,26.6249746504,26.1540541981,18.3553907645,24.8849430626] mm")
call updateGlobalVariable(oProject, "y_pts", "[0,0.75,0.771201413428,3.49095919837,6.97197534562,0,0,12.4233141649,7.3245789764,0,0,0,0,0,6.8819128783,8.49796734531,5.96402799144,7.03414591642] mm")
call updateGlobalVariable(oProject, "e_angles", "[0,0,0,-9.28748835642,0,15,0,-13.481382813,0,0,1.51861718696,6.71677945498,0,0,14.4923878983,3.50761210165,0,-18,0,0] deg")
