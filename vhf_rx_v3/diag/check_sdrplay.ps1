Get-Service -Name 'SDRplay*' | Select-Object Name, Status | Format-Table -AutoSize
Write-Output '---- PnP devices ----'
Get-PnpDevice -PresentOnly | Where-Object { $_.FriendlyName -like '*SDRplay*' -or $_.FriendlyName -like '*RSP*' } | Select-Object FriendlyName, Class, Status | Format-Table -AutoSize
Write-Output '---- all USB present (names) ----'
Get-PnpDevice -PresentOnly -Class USB | Select-Object -ExpandProperty FriendlyName
