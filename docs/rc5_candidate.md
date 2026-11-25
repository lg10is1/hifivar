# HiFiVar 0.1.0rc5 release identity / 发行身份

[简体中文](zh_CN/rc5_candidate.md)

[English installation](installation.md) · [中文安装](zh_CN/installation.md)

Use the [release assets](https://github.com/lg10is1/hifivar/releases/tag/v0.1.0-rc5) together with their checksums. Previous same-version development builds, including the 15 October payload, have different manifests and validation scopes. The public release source includes the indexed-BAM launcher, optional raw preparation, explicit source plans, support selection and attempt reporting. [Capabilities and evidence](capabilities.md)

公开使用以 GitHub 固定 tag 和校验表确认身份。早期同名开发包有各自的验证记录。当前功能范围见[中文说明](zh_CN/capabilities.md)，各研究脚本和论文资源按原身份追溯。

## Validated development payload / 已验收开发包

The **15 October 2026** development payload includes the bounded calling launcher,
raw-input preparation, offline review CLI, native preflight and opt-in
within-sample read-SV integration. Its Linux regression, package installation
and bounded real-tool gates have passed. The following hashes identify that
payload **before this documentation revision**, not every package named RC5:

| Artifact | Historical payload SHA256 |
| --- | --- |
| `hifivar-0.1.0rc5-py3-none-any.whl` | `3516cf88772eb8bbff1b12639d55cb732820d4d80f981fa23f72ce86eb14d4ab` |
| `hifivar-0.1.0rc5.tar.gz` | `70038ffeac144e816920bfbd1eb35537c0cef33deb78e7b59bbe829de816c9d5` |

Any build after a source or documentation change must carry a new delivery
manifest and measured artifact hashes. Do not copy these historical hashes
into a later `SHA256SUMS`, replace older assets or infer publication from a
successful local installation.

2026 年 10 月 15 日指定开发包已通过相应 Linux 验收，包含上述新增入口；最初的安全
修复包和中间开发包不因此变更。上表是**本轮文档修订之前的历史包身份**。源码或文档
改变后，交付清单和包哈希必须重新生成，不能沿用上表冒充字节一致。
