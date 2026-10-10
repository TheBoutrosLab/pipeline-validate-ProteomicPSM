import java.nio.charset.StandardCharsets
import java.nio.file.Files

/**
 * Expands environment variables in CSV inputs used by NFTest configurations.
 */
class NFTestCsvTools {
    /**
     * Expand $NAME and ${NAME} references using the environment inherited by
     * Nextflow. Undefined variables fail instead of leaving unusable paths.
     */
    String expandEnvironmentVariables(String content) {
        if (content == null) {
            return null
        }

        def variablePattern = ~/\$(?:\{([A-Za-z_][A-Za-z0-9_]*)\}|([A-Za-z_][A-Za-z0-9_]*))/

        return content.replaceAll(variablePattern) { match, bracedName, bareName ->
            def variableName = bracedName ?: bareName
            def variableValue = System.getenv(variableName)

            if (variableValue == null) {
                throw new IllegalArgumentException(
                    "Environment variable ${variableName} referenced by test input is not defined"
                )
            }

            return variableValue
        }
    }

    /**
     * Expand environment variables in a CSV and write the result beneath the
     * Nextflow work directory. NFTest removes this with the test work tree.
     */
    String prepareInputCsv(Object sourceCsv, String testName, Object launchDir) {
        def sourceFile = new File(sourceCsv.toString())
        if (!sourceFile.isFile()) {
            throw new IllegalArgumentException("Test input CSV does not exist: ${sourceFile}")
        }

        def workRoot = new File(System.getenv('NXF_WORK') ?: "${launchDir}/work")
        def inputDirectory = new File(workRoot, 'nftest-input')
        Files.createDirectories(inputDirectory.toPath())

        def safeTestName = testName.replaceAll(/[^A-Za-z0-9_.-]/, '-')
        def expandedCsv = File.createTempFile("${safeTestName}-", '.csv', inputDirectory)
        def sourceContent = Files.readString(sourceFile.toPath(), StandardCharsets.UTF_8)
        def expandedContent = expandEnvironmentVariables(sourceContent)

        Files.writeString(expandedCsv.toPath(), expandedContent, StandardCharsets.UTF_8)
        expandedCsv.deleteOnExit()

        return expandedCsv.absolutePath
    }
}

return new NFTestCsvTools()
