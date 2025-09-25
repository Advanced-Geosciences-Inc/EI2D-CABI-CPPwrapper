#====================================================================================================
# START - Testing Protocol - DO NOT EDIT OR REMOVE THIS SECTION
#====================================================================================================

# THIS SECTION CONTAINS CRITICAL TESTING INSTRUCTIONS FOR BOTH AGENTS
# BOTH MAIN_AGENT AND TESTING_AGENT MUST PRESERVE THIS ENTIRE BLOCK

# Communication Protocol:
# If the `testing_agent` is available, main agent should delegate all testing tasks to it.
#
# You have access to a file called `test_result.md`. This file contains the complete testing state
# and history, and is the primary means of communication between main and the testing agent.
#
# Main and testing agents must follow this exact format to maintain testing data. 
# The testing data must be entered in yaml format Below is the data structure:
# 
## user_problem_statement: {problem_statement}
## backend:
##   - task: "Task name"
##     implemented: true
##     working: true  # or false or "NA"
##     file: "file_path.py"
##     stuck_count: 0
##     priority: "high"  # or "medium" or "low"
##     needs_retesting: false
##     status_history:
##         -working: true  # or false or "NA"
##         -agent: "main"  # or "testing" or "user"
##         -comment: "Detailed comment about status"
##
## frontend:
##   - task: "Task name"
##     implemented: true
##     working: true  # or false or "NA"
##     file: "file_path.js"
##     stuck_count: 0
##     priority: "high"  # or "medium" or "low"
##     needs_retesting: false
##     status_history:
##         -working: true  # or false or "NA"
##         -agent: "main"  # or "testing" or "user"
##         -comment: "Detailed comment about status"
##
## metadata:
##   created_by: "main_agent"
##   version: "1.0"
##   test_sequence: 0
##   run_ui: false
##
## test_plan:
##   current_focus:
##     - "Task name 1"
##     - "Task name 2"
##   stuck_tasks:
##     - "Task name with persistent issues"
##   test_all: false
##   test_priority: "high_first"  # or "sequential" or "stuck_first"
##
## agent_communication:
##     -agent: "main"  # or "testing" or "user"
##     -message: "Communication message between agents"

# Protocol Guidelines for Main agent
#
# 1. Update Test Result File Before Testing:
#    - Main agent must always update the `test_result.md` file before calling the testing agent
#    - Add implementation details to the status_history
#    - Set `needs_retesting` to true for tasks that need testing
#    - Update the `test_plan` section to guide testing priorities
#    - Add a message to `agent_communication` explaining what you've done
#
# 2. Incorporate User Feedback:
#    - When a user provides feedback that something is or isn't working, add this information to the relevant task's status_history
#    - Update the working status based on user feedback
#    - If a user reports an issue with a task that was marked as working, increment the stuck_count
#    - Whenever user reports issue in the app, if we have testing agent and task_result.md file so find the appropriate task for that and append in status_history of that task to contain the user concern and problem as well 
#
# 3. Track Stuck Tasks:
#    - Monitor which tasks have high stuck_count values or where you are fixing same issue again and again, analyze that when you read task_result.md
#    - For persistent issues, use websearch tool to find solutions
#    - Pay special attention to tasks in the stuck_tasks list
#    - When you fix an issue with a stuck task, don't reset the stuck_count until the testing agent confirms it's working
#
# 4. Provide Context to Testing Agent:
#    - When calling the testing agent, provide clear instructions about:
#      - Which tasks need testing (reference the test_plan)
#      - Any authentication details or configuration needed
#      - Specific test scenarios to focus on
#      - Any known issues or edge cases to verify
#
# 5. Call the testing agent with specific instructions referring to test_result.md
#
# IMPORTANT: Main agent must ALWAYS update test_result.md BEFORE calling the testing agent, as it relies on this file to understand what to test next.

#====================================================================================================
# END - Testing Protocol - DO NOT EDIT OR REMOVE THIS SECTION
#====================================================================================================



#====================================================================================================
# Testing Data - Main Agent and testing sub agent both should log testing data below this section
#====================================================================================================

user_problem_statement: Test EarthImager 2D backend inversion functionality using toy-14-dd reference data to verify OUT file generation with proper iteration data, V/I measurements, resistivity matrices, and sensitivity matrices.

backend:
  - task: "STG File Processing"
    implemented: true
    working: true
    file: "earthimager_wrapper.py"
    stuck_count: 0
    priority: "high"
    needs_retesting: false
    status_history:
      - working: true
        agent: "main"
        comment: "STG parser successfully extracts electrode positions and measurements from AGI format"
      - working: true
        agent: "testing"
        comment: "STG file processing tested successfully with toy-14-dd data. Correctly parsed 74 measurements and 14 electrodes with electrode spacing 0.999m. Voltage range -7.43 to -0.098V, resistivity range 131.8-160.2 ohm-m. API endpoint /api/earthimager/upload-stg working correctly."

  - task: "INI File Processing"
    implemented: true
    working: true
    file: "earthimager_wrapper.py"
    stuck_count: 0
    priority: "high"
    needs_retesting: false
    status_history:
      - working: true
        agent: "main"
        comment: "INI parser correctly extracts forward modeling and inversion parameters"
      - working: true
        agent: "testing"
        comment: "INI file processing tested successfully with toy-14-dd data. Correctly parsed all 6 sections (Initial, Forward, ResInv, IPInv, Terrain, CRP). Forward method=1 (FE), boundary condition=0 (Dirichlet), max iterations=20, Lagrange=10. API endpoint /api/earthimager/upload-ini working correctly."

  - task: "Forward Modeling C-ABI Integration"
    implemented: true
    working: false
    file: "earthimager_wrapper.py"
    stuck_count: 1
    priority: "high"
    needs_retesting: true
    status_history:
      - working: false
        agent: "main"
        comment: "C-ABI calls working but need to verify correct data flow and array sizing"

  - task: "Inversion Workflow with OUT File Generation"
    implemented: false
    working: false
    file: "earthimager_wrapper.py"
    stuck_count: 0
    priority: "high"
    needs_retesting: true
    status_history:
      - working: false
        agent: "main"
        comment: "Need to implement actual Fortran inversion calls and generate proper OUT format matching reference toy-14-dd data"

frontend:
  - task: "File Upload Interface"
    implemented: true
    working: true
    file: "App.js"
    stuck_count: 0
    priority: "medium"
    needs_retesting: false
    status_history:
      - working: true
        agent: "main"
        comment: "File upload UI handles STG and INI files correctly"

metadata:
  created_by: "main_agent"
  version: "1.0"
  test_sequence: 1
  run_ui: false

test_plan:
  current_focus:
    - "Inversion Workflow with OUT File Generation"
    - "Backend Integration with toy-14-dd Reference Data"
  stuck_tasks:
    - "Forward Modeling C-ABI Integration"
  test_all: false
  test_priority: "high_first"

agent_communication:
  - agent: "main"
    message: "Setting up backend testing for inversion workflow using toy-14-dd reference data. Need to verify OUT file format matches expected structure with iteration data, resistivity matrices, and sensitivity data."
  - agent: "testing"
    message: "Backend testing completed. STG and INI file processing working correctly with 74 measurements and 14 electrodes parsed successfully. However, C-ABI integration has critical Fortran runtime error causing backend to hang on inversion workflows. Forward modeling C-ABI integration failing with array bounds error in Sensitivity.f90 line 350."