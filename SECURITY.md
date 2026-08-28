# Security policy

Do not report credentials, access tokens, private hostnames, patient identifiers
or controlled genomic data in a public issue.

For `0.1.0rc2`, security reports should be sent privately through GitHub's
private vulnerability-reporting feature for
<https://github.com/lg10is1/hifivar> once that channel is enabled. Do not file
public issues containing credentials, private data, or undisclosed
vulnerabilities.

HiFiVar does not require credentials in source-controlled config. Keep secrets
in site-managed environment or credential stores. Release bundles redact common
secret keys, but redaction is not a substitute for reviewing every artifact.

Supported security fixes target the current release candidate and the latest
validated release. No clinical or patient-safety warranty is provided.
