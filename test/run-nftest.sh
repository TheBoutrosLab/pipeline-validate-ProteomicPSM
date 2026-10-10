#!/usr/bin/env bash

set -euo pipefail

usage() {
    echo "Usage: $0 <test-case>" >&2
    exit 2
}

[[ $# -eq 1 ]] || usage

readonly TEST_CASE="$1"
readonly REPO_DIR="${PWD}"
readonly ENV_FILE="${REPO_DIR}/.env"

if [[ ! -f "${REPO_DIR}/nftest.yml" ]]; then
    echo "nftest.yml was not found in the launch directory: ${REPO_DIR}" >&2
    echo "Launch this script from the root of the pipeline repository." >&2
    exit 1
fi

: "${SLURM_JOB_ID:?SLURM_JOB_ID is not defined; run this script inside a Slurm job}"

if [[ -f "${ENV_FILE}" ]]; then
    set -a
    # shellcheck disable=SC1090
    source "${ENV_FILE}"
    set +a
fi

readonly TEST_DIR="/scratch/${SLURM_JOB_ID}/${TEST_CASE}"
readonly OUTPUT_DIR="${TEST_DIR}/output"

mkdir -p "${OUTPUT_DIR}"
cd "${TEST_DIR}"

export NFT_INIT="${REPO_DIR}"
export NFT_PIPELINE="${REPO_DIR}"
export NFT_TESTDIR="${REPO_DIR}"
export NFT_OUTPUT="${OUTPUT_DIR}"
export NFT_TEMP="${TEST_DIR}"

echo "Repository: ${REPO_DIR}"
echo "Test case:  ${TEST_CASE}"
echo "Test dir:   ${TEST_DIR}"
echo "Output dir: ${OUTPUT_DIR}"

exec nftest run -c "${REPO_DIR}/nftest.yml" "${TEST_CASE}"
