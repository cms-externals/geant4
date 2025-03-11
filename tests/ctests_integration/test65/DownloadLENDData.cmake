# - Script to download, unpack and move LEND data for test65
#
# Usage: cmake [-DALTERNATE_URL=<url>] -P DownloadLENDData.cmake
#
# When run, the script will check for a prexisting LEND data tarball
# under <CWD>/Data. If the file does not exist, or exists but has the
# wrong hash, download from LLNL is attempted. If the download fails,
# and ALTERNATE_URL has been supplied on the command line, download from
# ALTERNATE_URL/<filename> is attempted. In both cases, if download
# fails for any reason, a FATAL_ERROR is emitted.
#
# After successful download, the tarball is unpacked and its root
# directory relocated to <CWD>/Data/g4lend.
#
#
set(LENDDATA_FILENAME "LEND_GND1.3_ENDF.BVII.1.tar.gz")
set(LENDDATA_TAR_ROOTDIR "LEND_GND1.3_ENDF.BVII.1")
set(LENDDATA_HASH_MD5 "3e28151dc4c4647af3ae37d0385fc443")
set(LENDDATA_WORKING_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Data")
set(LENDDATA_LOCAL_FILENAME "${LENDDATA_WORKING_DIR}/${LENDDATA_FILENAME}")
set(LENDDATA_LOCAL_ROOTDIR "${LENDDATA_WORKING_DIR}/g4lend")

# Main and alternate URLs
set(LLNL_URL "ftp://gdo-nuclear.ucllnl.org/LEND_GND1.3/${LENDDATA_FILENAME}")

# If set, via command line arg
if(ALTERNATE_URL)
  string(REPLACE "\"" "" ALTERNATE_URL "${ALTERNATE_URL}" "/" "${LENDDATA_FILENAME}")
endif()

# Check for needed file, downloading if needed.
set(LENDDATA_NEEDS_DOWNLOAD TRUE)

if(EXISTS "${LENDDATA_LOCAL_FILENAME}")
  message(STATUS "test65-data: Checking MD5 hash of pre-existing file: ${LENDDATA_LOCAL_FILENAME}")
  file(MD5 "${LENDDATA_LOCAL_FILENAME}" __lenddata_existing_md5)
  if(__lenddata_existing_md5 STREQUAL ${LENDDATA_HASH_MD5})
    set(LENDDATA_NEEDS_DOWNLOAD FALSE)
    message(STATUS "MD5 hash ok: ${LENDDATA_LOCAL_FILENAME} = ${__lenddata_existing_md5}")
  endif()
endif()

if(LENDDATA_NEEDS_DOWNLOAD)
  message(STATUS "test65-data: No existing file, attempting download: ${LLNL_URL}")
  file(DOWNLOAD "${LLNL_URL}" "${LENDDATA_LOCAL_FILENAME}"
    SHOW_PROGRESS
    INACTIVITY_TIMEOUT 1200
    TIMEOUT 3000
    STATUS DownloadStatus
    )

  list(GET DownloadStatus 0 DownloadReturnStatus)
  if(DownloadReturnStatus AND ALTERNATE_URL)
    message(STATUS "test65-data: Download from ${LLNL_URL} failed, trying ${ALTERNATE_URL}")
    file(DOWNLOAD "${ALTERNATE_URL}" "${LENDDATA_LOCAL_FILENAME}"
      SHOW_PROGRESS
      INACTIVITY_TIMEOUT 600
      TIMEOUT 1500
      STATUS DownloadStatus
      )
  endif()

  list(GET DownloadStatus 0 DownloadReturnStatus)
  if(DownloadReturnStatus)
    message(FATAL_ERROR "test65-data: download FAILED: ${DownloadReturnStatus}, ${DownloadStringReturnStatus}")
  else()
    message(STATUS "test65-data: download OK")
  endif()

  message(STATUS "test65-data: Checking MD5 hash of download file: ${LENDDATA_LOCAL_FILENAME}")
  file(MD5 "${LENDDATA_LOCAL_FILENAME}" __lenddata_existing_md5)
  if(NOT __lenddata_existing_md5 STREQUAL ${LENDDATA_HASH_MD5})
    message(FATAL_ERROR "MD5 hash error: md5sum of download file = ${__lenddata_existing_md5}, expected ${LENDDATA_HASH_MD5}")
  endif()
endif()

# Remove any existing "g4lend" directory
if(EXISTS "${LENDDATA_LOCAL_ROOTDIR}")
  message(STATUS "Cleaning up existing unpack dir")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${LENDDATA_LOCAL_ROOTDIR}"
    WORKING_DIRECTORY "${LENDDATA_WORKING_DIR}"
    )
endif()

# Unpack downloaded tarball
message(STATUS "Attempting unpack: ${LENDDATA_LOCAL_FILENAME}")
execute_process(
  COMMAND ${CMAKE_COMMAND} -E tar xfz "${LENDDATA_LOCAL_FILENAME}"
  WORKING_DIRECTORY ${LENDDATA_WORKING_DIR}
  OUTPUT_QUIET
  RESULT_VARIABLE __lenddata_untar_result
  )

if(__lenddata_untar_result)
  message(FATAL_ERROR "test65-data: failed to untar file : ${LENDDATA_LOCAL_FILENAME}")
else()
  message(STATUS "test65-data: untarred '${LENDDATA_LOCAL_FILENAME}' OK")
endif()

# Rename versioned directory to "g4lend" for consistency
message(STATUS "Attempting rename of : ${LENDDATA_WORKING_DIR}/${LENDDATA_TAR_ROOTDIR} -> ${LENDDATA_LOCAL_ROOTDIR}")
execute_process(
  COMMAND ${CMAKE_COMMAND} -E rename "${LENDDATA_WORKING_DIR}/${LENDDATA_TAR_ROOTDIR}" "${LENDDATA_LOCAL_ROOTDIR}"
  WORKING_DIRECTORY ${LENDDATA_WORKING_DIR}
  OUTPUT_QUIET
  RESULT_VARIABLE __lenddata_rename_result
  )

if(__lenddata_rename_result)
  message(FATAL_ERROR "test65-data: failed to rename : ${LENDDATA_WORKING_DIR}/${LENDDATA_TAR_ROOTDIR} -> ${LENDDATA_LOCAL_ROOTDIR}")
else()
  message(STATUS "test65-data: renamed  OK")
endif()



