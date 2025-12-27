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


call updateGlobalVariable(oProject, "x_pts", "[0,80.9938269253,81.9937507145,82.0355140669,110.344043268,114,135,134.710954637,80.8265727823,81,82.0505119271,110.568475695,64.8749122227,71.7920677438,77.8082533063,72.3936863001,66.285684265,71.7002512712,74.7269399062,55,80.25,74.141332484,50.8133742881] mm")
call updateGlobalVariable(oProject, "y_pts", "[0,1,1.01234567901,1.56873811671,3.42417714774,0,0,8.82942244607,5.29765346764,0,0,0,0.5,0.5,17.0293339842,2.15293339842,4.37606433004,19.2524649158,27.5682235962,0,0,30.7103454473,21.0475887801] mm")
call updateGlobalVariable(oProject, "e_angles", "[0,0,0,-93.75,0,3.75,0,-3.04262672745,0,0,0.707373272549,1.09551479932,0,0,0,0,0,7.90474377685,0,0,0,0,22.5,0,-22.5] deg")
