# Presence of this file tells shiny::loadSupport() NOT to auto-source every
# R/*.R file when runApp() starts. CanProTarget loads helpers explicitly in
# app.R. Without this, mcp_worker.R (an MCP stdin server) would also be sourced
# at app boot, loading DepMap twice and briefly entering a worker main loop.
#
# See: shiny:::loadSupport() / shiny docs on the R/ support directory.
