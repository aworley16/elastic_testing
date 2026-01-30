#include <pmix.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    pmix_status_t rc;
    pmix_proc_t myproc;
    pmix_info_t *info = NULL;
    size_t ninfo = 0;
    
    pmix_info_t *results = NULL;
    size_t nresults = 0;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <target_namespace>\n", argv[0]);
        return 1;
    }

    const char *target_nspace = argv[1];

    /* Initialize PMIx */
    rc = PMIx_Init(&myproc, NULL, 0);
    if (rc != PMIX_SUCCESS) {
        fprintf(stderr, "PMIx_Init failed: %s\n", PMIx_Error_string(rc));
        return 1;
    }

    printf("PMIx client initialized as %s:%d\n",
           myproc.nspace, myproc.rank);

    /* Ask PMIx for job size of the target namespace */
    PMIX_QUERY_CONSTRUCT(&query);
    query.keys = (char **)malloc(2 * sizeof(char *));
    query.keys[0] = strdup(PMIX_JOB_SIZE);
    query.keys[1] = NULL;

    PMIX_INFO_CREATE(query.qualifiers, 1);
    PMIX_INFO_LOAD(&query.qualifiers[0],
                   PMIX_NSPACE, target_nspace, PMIX_STRING);
    query.nqual = 1;

    rc = PMIx_Query_info(&query, 1, &results, &nresults);
    if (rc != PMIX_SUCCESS) {
        fprintf(stderr, "PMIx_Query_info failed: %s\n",
                PMIx_Error_string(rc));
        goto cleanup;
    }

    for (size_t i = 0; i < nresults; i++) {
        if (strcmp(results[i].key, PMIX_JOB_SIZE) == 0) {
            uint32_t job_size = results[i].value.data.uint32;
            printf("Found MPI job %s with %u ranks\n",
                   target_nspace, job_size);

            printf("Processes in job:\n");
            for (uint32_t r = 0; r < job_size; r++) {
                printf("  %s:%u\n", target_nspace, r);
            }
        }
    }

cleanup:
    if (results) {
        PMIX_INFO_FREE(results, nresults);
    }
    PMIX_QUERY_DESTRUCT(&query);

    PMIx_Finalize(NULL, 0);
    return 0;
}
