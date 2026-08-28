//
//  Functions used by Geant4 pipeline jobs
//
//  Gunter Folger, April 2021
//

def getLABEL() {
  def LABEL = 'docker-host-noafs'
  if (params.OS =~ '^g4-' ) LABEL=params.OS
  if (params.OS =~ '^mac' ) LABEL=params.OS
  if (params.OS =~ 'arm64$' ) LABEL='arm64-docker'
  if (params.OS =~ 'windows10')  LABEL='g4-win10-' + params.COMPILER
  if (params.OS =~ 'windows11')  LABEL='g4-win11-' + params.COMPILER
  return LABEL
}

//- for g4- and mac nodes, set a second label used DoBuild() when running on hardware
def getLABELSelector() {
  def value = 'INVALID'
  if (params.OS =~ '^g4-' ) value='g4-node'
  if (params.OS =~ '^mac' ) value='baremacel'
  return value
}

def GetDockerImage() {
  def IMAGE="sft/docker/" + "${OS}"
  if ("${OS}" =~ '^centos8' ) IMAGE="${IMAGE}" + "-g4"
  return "gitlab-registry.cern.ch/" + "${IMAGE}"
}


def GetBuildOptionDir() {
  def dir=params.BuildOptions.replaceAll(' ','').replaceAll(',','_')
  if (params.BuildOptions.equals('none')) dir=''
  return dir
}

def GetBuildOptionText() {
  def txt='-' + GetBuildOptionDir()
  if (txt.equals('-')) txt=''
  return txt
}

def GetToday() {
  def txt="${new Date().format('yyyy-MM-dd')}"
  return txt
}

def doCheckout(needVecGeom) {
  checkout changelog: false,
    poll: false,
    scm: [$class: 'GitSCM', branches: [[name: '${SourceBranch}']],
          extensions: [[$class: 'CheckoutOption', timeout: 15], [$class: 'RelativeTargetDirectory', relativeTargetDir: 'geant4'],
                       [$class: 'CloneOption', honorRefspec: true, noTags: false, reference: '', shallow: false, timeout: 20]],
          userRemoteConfigs: [[credentialsId: 'ade6d0b4-c04e-401f-b04b-1e6ef81c5ed7',
                        refspec: "+refs/heads/master:refs/remotes/origin/master +refs/${SourceBranch}/latest:refs/${SourceBranch}",
          url: "${SourceRepoHttpUrl}"]]]

  checkout changelog: false,
    poll: false,
    scm: [$class: 'GitSCM', branches: [[name: '${VerificationSourceBranch}']],
          extensions: [[$class: 'CheckoutOption', timeout: 15], [$class: 'RelativeTargetDirectory', relativeTargetDir: 'geant4/verification'],
                       [$class: 'CloneOption', honorRefspec: true, noTags: false, reference: '', shallow: false, timeout: 14]],
          userRemoteConfigs: [[credentialsId: 'ade6d0b4-c04e-401f-b04b-1e6ef81c5ed7',
                          refspec: "+refs/heads/master:refs/remotes/origin/master +refs/${VerificationSourceBranch}/latest:refs/${VerificationSourceBranch}",
                          url: "${VerificationRepoHttpUrl}"]]]

  checkout changelog: false,
    poll: false,
    scm: [$class: 'GitSCM', branches: [[name: '${BenchmarksSourceBranch}']],
          extensions: [[$class: 'CheckoutOption', timeout: 15], [$class: 'RelativeTargetDirectory', relativeTargetDir: 'geant4/benchmarks'],
                       [$class: 'CloneOption', honorRefspec: true, noTags: false, reference: '', shallow: false, timeout: 14]],
          userRemoteConfigs: [[credentialsId: 'ade6d0b4-c04e-401f-b04b-1e6ef81c5ed7',
                          refspec: "+refs/heads/master:refs/remotes/origin/master +refs/${BenchmarksSourceBranch}/latest:refs/${BenchmarksSourceBranch}",
                          url: "${BenchmarksRepoHttpUrl}"]]]

  if ( needVecGeom ) {
    checkout([$class: 'GitSCM', branches: [[name: '*/master']],
              extensions: [[$class: 'RelativeTargetDirectory', relativeTargetDir: 'geant4/VecGeom']],
              userRemoteConfigs: [[url: 'ssh://git@gitlab.cern.ch:7999/VecGeom/VecGeom.git']]])
  }
}


def doDockerBuildTest() {
  sh 'touch controlfile'
  sh 'date'
  sh 'hostname'
  sh 'env | sort'
  sh 'date'
  withEnv(["build_script=${scriptToRun}"]) {
    sh 'scripts/tests/tools/ctest/docker-runbuild.sh'
  }
  sh 'date'
  sh 'test "$MODE" != "continuous" && scripts/tests/tools/ctest/copy_to_eos.sh || true'
  sh 'date'
}

def doHWBuildTest() {
  sh 'touch controlfile'
  sh 'date'
  sh 'hostname'
  sh 'env | sort'
  sh 'date'
  withEnv(["LABEL=${OS}"]) {
    sh 'scripts/tests/tools/ctest/jk-all'
  }
  sh 'date'
}

def doWinBuildTest() {
  bat 'echo %PATH%'
  bat 'echo BuildOptions=%BuildOptions%'
  bat 'python scripts/tests/tools/ctest/g4jk-win-run-pipe.py'
}

//-----------
def doBuild() {
  println "submit job for " + OS + " : " + COMPILER + " : " + THREAD + " : " + BuildOptions
  build wait: false,
    job: '/g4np-build',
    parameters: [
      string(name: 'OS', value: "${OS}"),
      string(name: 'COMPILER', value: "${COMPILER}"),
      string(name: 'THREAD', value: "${THREAD}"),
      string(name: 'MODE', value: "$MODE"),
      string(name: 'BUILDTYPE', value: "${BUILDTYPE}"),
      string(name: 'EXTERNALS', value: "${EXTERNALS}"),
      string(name: 'BuildOptions', value: "${BuildOptions}"),
      string(name: 'ExtraCMakeOptions', value: "${ExtraCMakeOptions}"),
      string(name: 'SourceBranch', value: "${SourceBranch}"),
      string(name: 'VerificationSourceBranch', value: "${VerificationSourceBranch}"),
      string(name: 'BenchmarksSourceBranch', value: "${BenchmarksSourceBranch}"),
    ]
}

def doOptionalBuilds(OS, compiler, thread, boList){
  for (opt in boList) {
    println "submit job for " + OS + " : " + compiler + " : " + thread + " : " + opt
    build wait: false,
      job: '/g4np-build',
      parameters: [
        string(name: 'OS', value: OS),
        string(name: 'COMPILER', value: compiler),
        string(name: 'THREAD', value: thread),
        string(name: 'MODE', value: "$MODE"),
        string(name: 'BUILDTYPE', value: "${BUILDTYPE}"),
        string(name: 'EXTERNALS', value: "${EXTERNALS}"),
        string(name: 'BuildOptions', value: opt),
        string(name: 'ExtraCMakeOptions', value: "${ExtraCMakeOptions}"),
        string(name: 'SourceBranch', value: "${SourceBranch}"),
        string(name: 'VerificationSourceBranch', value: "${VerificationSourceBranch}"),
        string(name: 'BenchmarksSourceBranch', value: "${BenchmarksSourceBranch}"),
      ]
  }
}

def doContinuousBuild(OS, compiler, thread, buildtype, boList) {
  for (opt in boList) {
    println "submit job for " + OS + " : " + compiler + " : " + thread + " : " + opt
    build wait: true,
      job: '/g4p-CI-build',
      parameters: [
        string(name: 'OS', value: OS),
        string(name: 'COMPILER', value: compiler),
        string(name: 'THREAD', value: thread),
        string(name: 'BUILDTYPE', value: buildtype),
        string(name: 'EXTERNALS', value: "${EXTERNALS}"),
        string(name: 'BuildOptions', value: opt),
        string(name: 'TargetRepoUrl', value: "${TargetRepoUrl}"),
        string(name: 'TargetBranch', value: "${TargetBranch}"),
        string(name: 'MergeRequestId', value: "${MergeRequestId}"),
        string(name: 'MergeRequestLastCommit', value: "${MergeRequestLastCommit}"),
        string(name: 'gitlabUserName', value: "${gitlabUserName}"),
      ]
  }
}

def abortPreviousCI() {
  def jobName = env.JOB_NAME
  def buildNumber = env.BUILD_NUMBER.toInteger()
  def currentJob = Jenkins.instance.getItemByFullName(jobName)

  for (def build : currentJob.builds) {
    def exec = build.getExecutor()
    println "abortPreviousCI debug: build.displayName= " + build.displayName
    if (build.displayName =~ " MR-${MergeRequestId}" ) {
      if (build.isBuilding() && build.number.toInteger() != buildNumber && exec != null) {
        exec.interrupt(
          Result.ABORTED,
          new CauseOfInterruption.UserInterruption("Aborted by #${currentBuild.number}")
        )
        println "===>>> Aborted ${build.number} for MR ${MergeRequestId}"
      }
   }
  }
}

def CheckForCvmfs() {
  if (params.OS =~ '^mac' ) {
     sh 'scripts/tests/tools/ctest/CheckForCvmfs.sh'
  }
}

return this
