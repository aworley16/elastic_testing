#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    MPI_Session session;
    MPI_Info info;
    int rc;

    if (argc < 3) {
        fprintf(stderr,
                "Usage: %s <target_pmix_namespace> <job_size>\n",
                argv[0]);
        return 1;
    }

    const char *target_nspace = argv[1];
    int job_size = atoi(argv[2]);

    /* ----------------------------------------
     * 1. Initialize MPI Session
     * ---------------------------------------- */
    MPI_Info_create(&info);

    rc = MPI_Session_init(info, MPI_ERRORS_RETURN, &session);
    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Session_init failed\n");
        return 1;
    }

    printf("MPI Session initialized\n");

    /* ----------------------------------------
     * 2. Build pset specification strings
     *    Format: mpi://<nspace>:<rank>
     * ---------------------------------------- */
    char **pset_specs = malloc(job_size * sizeof(char *));
    if (!pset_specs) {
        fprintf(stderr, "Out of memory\n");
        MPI_Session_finalize(&session);
        return 1;
    }

    for (int i = 0; i < job_size; i++) {
        /* room for "mpi://" + nspace + ":" + rank */
        pset_specs[i] = malloc(256);
        snprintf(pset_specs[i], 256,
                 "mpi://%s:%d",
                 target_nspace, i);
    }

    /* ----------------------------------------
     * 3. Create the dynamic pset
     * ---------------------------------------- */
    const char *new_pset_name = "dynamic_worker_pset";

    rc = MPI_Session_pset_create(session,
                                 new_pset_name,
                                 MPI_INFO_NULL,
                                 (const char **)pset_specs,
                                 job_size);

    if (rc != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Session_pset_create failed\n");
    } else {
        printf("Created new pset '%s' with %d processes\n",
               new_pset_name, job_size);
    }

    /* ----------------------------------------
     * 4. (Optional) Query the pset size
     * ---------------------------------------- */
    int pset_size;
    rc = MPI_Session_get_pset_info(session,
                                   new_pset_name,
                                   &info);

    if (rc == MPI_SUCCESS) {
        char size_str[32];
        int flag;
        MPI_Info_get(info, "mpi_pset_size",
                     sizeof(size_str),
                     size_str, &flag);
        if (flag) {
            pset_size = atoi(size_str);
            printf("Verified pset size: %d\n", pset_size);
        }
        MPI_Info_free(&info);
    }

    /* ----------------------------------------
     * 5. Cleanup
     * ---------------------------------------- */
    for (int i = 0; i < job_size; i++) {
        free(pset_specs[i]);
    }
    free(pset_specs);

    MPI_Session_finalize(&session);
    return 0;
}
