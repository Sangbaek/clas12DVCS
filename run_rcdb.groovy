import run.Run

output = ""
new File("/work/clas12/sangbaek/inb_run_list_pass1").eachLine{ line ->
  int runnum = line.toInteger()
  def run = new Run(runnum)
  beam_current_request = run.rcdb.rcdb_dict["beam_current_request"]
  
  output = output + line + ","+ beam_current_request + "\n"
}
new File("/work/clas12/sangbaek/outb_run_list_pass1").eachLine{ line ->
  int runnum = line.toInteger()
  def run = new Run(runnum)
  beam_current_request = run.rcdb.rcdb_dict["beam_current_request"]
  
  output = output + line + ","+ beam_current_request + "\n"
}

//println(run.special_runs)
def myFile = new File("rcdb.csv")
myFile.write(output)
