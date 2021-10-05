#!/usr/bin/env python

import json
import sys
from globus_sdk import (NativeAppAuthClient,
                        RefreshTokenAuthorizer,
                        TransferClient,
                        TransferData)
from globus_sdk.exc import GlobusAPIError


CLIENT_ID = "a5f78945-0459-4339-94ba-3e68e00b011e"
TOKEN_FILE = "refresh-tokens.json"

NMDC_ENDPOINT_ID = "72dd396a-2242-11ec-a0a4-6b21ca6daf73"
KBASE_ENDPOINT_ID = "c3c0a65f-5827-4834-b6c9-388b0b19953a"


def load_tokens_from_file(filepath):
    """Load a set of saved tokens."""
    with open(filepath, "r") as f:
        tokens = json.load(f)

    return tokens


def save_tokens_to_file(filepath, tokens):
    """Save a set of tokens for later use."""
    with open(filepath, "w") as f:
        json.dump(tokens, f)


def update_tokens_file_on_refresh(token_response):
    """
    Callback function passed into the RefreshTokenAuthorizer.
    Will be invoked any time a new access token is fetched.
    """
    save_tokens_to_file(TOKEN_FILE, token_response.by_resource_server)


def transfer(flist):
    tokens = load_tokens_from_file(TOKEN_FILE)

    transfer_tokens = tokens["transfer.api.globus.org"]

    auth_client = NativeAppAuthClient(client_id=CLIENT_ID)

    authorizer = RefreshTokenAuthorizer(
        transfer_tokens["refresh_token"],
        auth_client,
        access_token=transfer_tokens["access_token"],
        expires_at=transfer_tokens["expires_at_seconds"],
        on_refresh=update_tokens_file_on_refresh,
    )

    transfer = TransferClient(authorizer=authorizer)

    try:
        transfer.endpoint_autoactivate(NMDC_ENDPOINT_ID)
        transfer.endpoint_autoactivate(KBASE_ENDPOINT_ID)
    except GlobusAPIError as ex:
        print(ex)
        if ex.http_status == 401:
            sys.exit(
                "Refresh token has expired. "
                "Please delete refresh-tokens.json and try again."
            )
        else:
            raise ex

    tdata = TransferData(transfer,
                         NMDC_ENDPOINT_ID,
                         KBASE_ENDPOINT_ID,
                         label="Automated Transfers",
                         sync_level="checksum")
    for srcfile in flist:
        fn = srcfile.split('/')[-1]
        tdata.add_item(srcfile, "/%s/%s" % ("nmdc", fn))
    transfer_result = transfer.submit_transfer(tdata)
    print("task_id =", transfer_result["task_id"])


if __name__ == "__main__":
    transfer(sys.argv[1:])
