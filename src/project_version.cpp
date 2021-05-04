
namespace MR {
  namespace App {
    extern const char* project_version;
    extern const char* project_build_date;
    void set_project_version () {
      project_version = "unknown";
      project_build_date = __DATE__;
    }
  }
}
